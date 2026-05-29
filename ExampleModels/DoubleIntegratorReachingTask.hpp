#pragma once

#include <array>
#include <cstddef>
#include <iostream>

#include "AugmentedLagrangian/StateAugmentedLagrangian.hpp"
#include "Constraint/LinearStateConstraint.hpp"
#include "Cost/QuadraticStateCost.hpp"
#include "Dynamics/LinearSystemDynamics.hpp"
#include "Initialization/DefaultInitializer.hpp"
#include "OptimalControl/OptimalControlProblem.hpp"
#include "Penalties/QuadraticPenalty.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"
#include "iLQR/iLQRDescriptor.hpp"

namespace double_integrator {
template <typename Scalar, int ArrayLength>
class DoubleIntegratorReachingTask {
 public:
  static_assert(ArrayLength > 1,
                "DoubleIntegratorReachingTask requires at least two nodes.");

  static constexpr int STATE_DIM = 2;
  static constexpr int INPUT_DIM = 1;
  static constexpr int GOAL_CONSTRAINT_DIM = STATE_DIM;
  static constexpr Scalar timeStep = 1e-2;
  static constexpr Scalar minRelCost = 1e-12;  // 避免过早终止。
  static constexpr Scalar constraintTolerance = 1e-3;
  using StateVector_t = Vector<Scalar, STATE_DIM>;
  using StateMatrix_t = Matrix<Scalar, STATE_DIM, STATE_DIM>;
  using InputVector_t = Vector<Scalar, INPUT_DIM>;
  using InputMatrix_t = Matrix<Scalar, INPUT_DIM, INPUT_DIM>;
  using StateInputMatrix_t = Matrix<Scalar, STATE_DIM, INPUT_DIM>;
  using TimeTrajectory_t = std::array<Scalar, ArrayLength>;
  using StateTrajectory_t = std::array<StateVector_t, ArrayLength>;
  using InputTrajectory_t = std::array<InputVector_t, ArrayLength>;
  using DefaultInitializer_t = DefaultInitializer<Scalar, STATE_DIM, INPUT_DIM>;
  using QuadraticStateInputCost_t =
      QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>;
  using LinearSystemDynamics_t =
      LinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM>;
  using GoalConstraint_t =
      LinearStateConstraint<Scalar, STATE_DIM, GOAL_CONSTRAINT_DIM>;
  using GoalPenalty_t = QuadraticPenalty<Scalar>;
  using GoalAugmentedLagrangian_t =
      StateAugmentedLagrangian<Scalar, STATE_DIM, GOAL_CONSTRAINT_DIM>;

  DoubleIntegratorReachingTask() = default;
  virtual ~DoubleIntegratorReachingTask() = default;

 protected:
  const Scalar tGoal = 1.0;
  const StateVector_t xInit = StateVector_t::Zero();
  const StateVector_t xGoal{Scalar(2.0), Scalar(0.0)};

 public:
  class TerminalGoalAugmentedLagrangian final {
   public:
    explicit TerminalGoalAugmentedLagrangian(const StateVector_t& target)
        : constraint_(-target, StateMatrix_t::Identity()),
          scalarPenalty_(
              typename GoalPenalty_t::Config{penaltyScale, stepSize}),
          augmentedLagrangian_(&constraint_, &scalarPenalty_) {}

    GoalAugmentedLagrangian_t* get() { return &augmentedLagrangian_; }
    const GoalAugmentedLagrangian_t* get() const {
      return &augmentedLagrangian_;
    }

    const GoalConstraint_t& constraint() const { return constraint_; }

   private:
    static constexpr Scalar penaltyScale = 10.0;
    static constexpr Scalar stepSize = 1.0;

    GoalConstraint_t constraint_;
    GoalPenalty_t scalarPenalty_;
    GoalAugmentedLagrangian_t augmentedLagrangian_;
  };

  /*
   * 该任务有意保持为单模式：
   *   min integral 0.5 * u'Ru dt。
   *   s.t. x_dot = [v, u], x(0) = xInit, x(tGoal) = xGoal。
   *
   * 将 TerminalGoalAugmentedLagrangian::get() 绑定到
   * problem.finalEqualityLagrangian。不使用中间模式相关的
   * 状态-输入约束。
   */
  void configureTargetTrajectory(TimeTrajectory_t& timeTrajectory,
                                 StateTrajectory_t& stateTrajectory,
                                 InputTrajectory_t& inputTrajectory) const {
    for (size_t i = 0; i < timeTrajectory.size(); ++i) {
      const Scalar alpha = static_cast<Scalar>(i) /
                           static_cast<Scalar>(timeTrajectory.size() - 1);
      timeTrajectory[i] = alpha * tGoal;
      stateTrajectory[i] = (Scalar(1) - alpha) * xInit + alpha * xGoal;
      inputTrajectory[i].setZero();
    }
  }

  /*
   * 打印轨迹。可在 MATLAB 中使用以下命令绘图：
   * subplot(2, 1, 1); plot(timeTrajectory, stateTrajectory);
   * xlabel("time [sec]"); legend("pos", "vel");
   * subplot(2, 1, 2); plot(timeTrajectory, inputTrajectory);
   * xlabel("time [sec]"); legend("force");
   */
  template <typename PrimalSolution_t>
  void printSolution(const PrimalSolution_t& primalSolution,
                     bool display) const {
    if (display) {
      std::cerr << "\n";
      // 时间
      std::cerr << "timeTrajectory = [";
      for (const auto& t : primalSolution.timeTrajectory_) {
        std::cerr << t << "; ";
      }
      std::cerr << "];\n";
      // 状态
      std::cerr << "stateTrajectory = [";
      for (const auto& x : primalSolution.stateTrajectory_) {
        std::cerr << x.transpose() << "; ";
      }
      std::cerr << "];\n";
      // 输入。
      std::cerr << "inputTrajectory = [";
      for (const auto& u : primalSolution.inputTrajectory_) {
        std::cerr << u.transpose() << "; ";
      }
      std::cerr << "];\n";
    }
  }
};

template <typename Scalar, std::size_t PredictLength>
using DoubleIntegratorReachingTaskForHorizon =
    DoubleIntegratorReachingTask<Scalar, static_cast<int>(PredictLength + 1)>;

template <typename Scalar, std::size_t PredictLength>
using DoubleIntegratorReachingConstraintConfig =
    ConstraintConfig<StateConstraintConfig<ConstraintLayout<>>,
                     StateInputConstraintConfig<ConstraintLayout<>>,
                     FinalStateConstraintConfig<ConstraintLayout<
                         ConstraintGroupLayout<ConstraintTerm<
                             DoubleIntegratorReachingTaskForHorizon<
                                 Scalar, PredictLength>::GOAL_CONSTRAINT_DIM>>,
                         ConstraintGroupLayout<>>>>;

template <typename Scalar, std::size_t PredictLength>
using DoubleIntegratorReachingProblem = OptimalControlProblem<
    Scalar,
    TranscriptionConfig<Dimensions<DoubleIntegratorReachingTaskForHorizon<
                                       Scalar, PredictLength>::STATE_DIM,
                                   DoubleIntegratorReachingTaskForHorizon<
                                       Scalar, PredictLength>::INPUT_DIM>,
                        Horizon<PredictLength>>,
    DoubleIntegratorReachingConstraintConfig<Scalar, PredictLength>>;

template <typename Scalar, std::size_t PredictLength>
inline DoubleIntegratorReachingProblem<Scalar, PredictLength>&
createDoubleIntegratorReachingProblem() {
  using Task_t = DoubleIntegratorReachingTaskForHorizon<Scalar, PredictLength>;
  using StateVector_t = typename Task_t::StateVector_t;
  using StateMatrix_t = typename Task_t::StateMatrix_t;
  using InputMatrix_t = typename Task_t::InputMatrix_t;
  using StateInputMatrix_t = typename Task_t::StateInputMatrix_t;
  using Dynamics_t = typename Task_t::LinearSystemDynamics_t;
  using Cost_t = typename Task_t::QuadraticStateInputCost_t;
  using TerminalGoal_t = typename Task_t::TerminalGoalAugmentedLagrangian;
  using Problem_t = DoubleIntegratorReachingProblem<Scalar, PredictLength>;

  static const StateMatrix_t A{{Scalar(0.0), Scalar(1.0)},
                               {Scalar(0.0), Scalar(0.0)}};
  static const StateInputMatrix_t B{Scalar(0.0), Scalar(1.0)};
  static Dynamics_t dynamics(A, B);

  static const StateMatrix_t Q = StateMatrix_t::Zero();
  static const InputMatrix_t R = Scalar(0.1) * InputMatrix_t::Identity();
  static Cost_t cost(Q, R, 0);

  static const StateVector_t xGoal{Scalar(2.0), Scalar(0.0)};
  static TerminalGoal_t terminalGoal(xGoal);
  static Problem_t problem;
  static bool initialized = false;

  if (!initialized) {
    problem.dynamicsPtr = &dynamics;
    problem.cost.add(cost);
    problem.finalEqualityLagrangian.template set<0>(terminalGoal.get());
    initialized = true;
  }

  return problem;
}
};  // namespace double_integrator
