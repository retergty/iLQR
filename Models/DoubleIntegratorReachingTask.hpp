#pragma once

#include <cstddef>
#include <iostream>
#include <memory>

#include "DefaultInitializer.hpp"
#include "LinearSystemDynamics.hpp"
#include "OptimalControlProblem.hpp"
#include "PrimalSolution.hpp"
#include "QuadraticPenalty.hpp"
#include "QuadraticStateCost.hpp"
#include "SmoothAbsolutePenalty.hpp"

namespace double_integrator {
template <typename Scalar, int ArrayLength> class DoubleIntegratorReachingTask {
public:
  static constexpr int STATE_DIM = 2;
  static constexpr int INPUT_DIM = 1;
  static constexpr Scalar timeStep = 1e-2;
  static constexpr Scalar minRelCost = 1e-12; // to avoid early termination
  static constexpr Scalar constraintTolerance = 1e-3;
  using StateVector_t = Vector<Scalar, STATE_DIM>;
  using StateMatrix_t = Matrix<Scalar, STATE_DIM, STATE_DIM>;
  using InputVector_t = Vector<Scalar, INPUT_DIM>;
  using InputMatrix_t = Matrix<Scalar, INPUT_DIM, INPUT_DIM>;
  using StateInputMatrix_t = Matrix<Scalar, STATE_DIM, INPUT_DIM>;
  using DefaultInitializer_t = DefaultInitializer<Scalar, STATE_DIM, INPUT_DIM>;
  using StateInputCost_t = StateInputCost<Scalar, STATE_DIM, INPUT_DIM>;
  using QuadraticStateInputCost_t =
      QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>;
  using LinearSystemDynamics_t =
      LinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM>;
  using SystemDynamicsBase_t = SystemDynamicsBase<Scalar, STATE_DIM, INPUT_DIM>;

  enum class PenaltyType {
    QuadraticPenalty,
    SmoothAbsolutePenalty,
  };

  DoubleIntegratorReachingTask() = default;
  virtual ~DoubleIntegratorReachingTask() = default;

protected:
  const Scalar tGoal = 1.0;
  const StateVector_t xInit = StateVector_t::Zero();
  const StateVector_t xGoal = (StateVector_t << 2.0, 0.0).finished();

  std::unique_ptr<DefaultInitializer_t> getInitializer() const {
    return std::make_unique<DefaultInitializer_t>();
  }

  std::shared_ptr<ReferenceManager> getReferenceManagerPtr() const {
    const ModeSchedule modeSchedule({tGoal}, {0, 1});
    const TargetTrajectories targetTrajectories({tGoal}, {xGoal},
                                                {vector_t::Zero(INPUT_DIM)});
    return std::make_shared<ReferenceManager>(targetTrajectories, modeSchedule);
  }

  std::unique_ptr<StateInputCost_t> getCostPtr() const {
    StateMatrix_t Q = StateMatrix_t::Zero();
    InputMatrix_t R = 0.1 * InputMatrix_t::Identity();
    return std::make_unique<QuadraticStateInputCost_t>(std::move(Q),
                                                       std::move(R));
  }

  std::unique_ptr<SystemDynamicsBase_t> getDynamicsPtr() const {
    StateMatrix_t A = (StateMatrix_t() << 0.0, 1.0, 0.0, 0.0).finished();
    StateInputMatrix_t B = (StateInputMatrix_t << 0.0, 1.0).finished();
    return std::make_unique<LinearSystemDynamics_t>(std::move(A), std::move(B));
  }

  std::unique_ptr<StateAugmentedLagrangian>
  getGoalReachingAugmentedLagrangian(const StateVector_t &xGoal,
                                     PenaltyType penaltyType) {
    constexpr Scalar scale = 10.0;
    constexpr Scalar stepSize = 1.0;

    auto goalReachingConstraintPtr = std::make_unique<LinearStateConstraint>(
        -xGoal, matrix_t::Identity(xGoal.size(), xGoal.size()));

    switch (penaltyType) {
    case PenaltyType::QuadraticPenalty: {
      using penalty_type = augmented::QuadraticPenalty;
      penalty_type::Config boundsConfig{scale, stepSize};
      return create(std::move(goalReachingConstraintPtr),
                    penalty_type::create(boundsConfig));
    }
    case PenaltyType::SmoothAbsolutePenalty: {
      using penalty_type = augmented::SmoothAbsolutePenalty;
      penalty_type::Config boundsConfig{scale, 0.1, stepSize};
      return create(std::move(goalReachingConstraintPtr),
                    penalty_type::create(boundsConfig));
    }
    default:
      return nullptr;
    }
  }

  /** This class enforces zero force at the second mode (mode = 1)*/
  class ZeroInputConstraint final : public StateInputConstraint {
  public:
    ZeroInputConstraint(const ReferenceManager &referenceManager)
        : StateInputConstraint(ConstraintOrder::Linear),
          referenceManagerPtr_(&referenceManager) {}

    ~ZeroInputConstraint() override = default;
    ZeroInputConstraint *clone() const override {
      return new ZeroInputConstraint(*this);
    }

    size_t getNumConstraints(Scalar time) const override { return 1; }

    /** Only active after the fist mode */
    bool isActive(Scalar time) const override {
      return referenceManagerPtr_->getModeSchedule().modeAtTime(time) > 0
                 ? true
                 : false;
    }

    vector_t getValue(Scalar t, const vector_t &x, const vector_t &u,
                      const PreComputation &) const override {
      return u;
    }

    VectorFunctionLinearApproximation
    getLinearApproximation(Scalar t, const vector_t &x, const vector_t &u,
                           const PreComputation &) const override {
      VectorFunctionLinearApproximation approx;
      approx.f = u;
      approx.dfdx.setZero(getNumConstraints(t), x.size());
      approx.dfdu.setIdentity(getNumConstraints(t), u.size());
      return approx;
    }

  private:
    ZeroInputConstraint(const ZeroInputConstraint &) = default;
    const ReferenceManager *referenceManagerPtr_;
  };

  /*
   * printout trajectory. Use the following commands for plotting in MATLAB:
   * subplot(2, 1, 1); plot(timeTrajectory, stateTrajectory);
   * xlabel("time [sec]"); legend("pos", "vel");
   * subplot(2, 1, 2); plot(timeTrajectory, inputTrajectory);
   * xlabel("time [sec]"); legend("force");
   */
  void printSolution(const PrimalSolution &primalSolution, bool display) const {
    if (display) {
      std::cerr << "\n";
      // time
      std::cerr << "timeTrajectory = [";
      for (const auto &t : primalSolution.timeTrajectory_) {
        std::cerr << t << "; ";
      }
      std::cerr << "];\n";
      // state
      std::cerr << "stateTrajectory = [";
      for (const auto &x : primalSolution.stateTrajectory_) {
        std::cerr << x.transpose() << "; ";
      }
      std::cerr << "];\n";
      // input
      std::cerr << "inputTrajectory = [";
      for (const auto &u : primalSolution.inputTrajectory_) {
        std::cerr << u.transpose() << "; ";
      }
      std::cerr << "];\n";
    }
  }
};
}; // namespace double_integrator
