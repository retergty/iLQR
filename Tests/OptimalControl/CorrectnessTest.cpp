#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

#include "DefaultInitializer.hpp"
#include "LinearSystemDynamics.hpp"
#include "MatrixEigenConversion.hpp"
#include "QpDiscreteTranscription.hpp"
#include "QpSolver.hpp"
#include "QuadraticPenalty.hpp"
#include "SensitivityIntegrator.hpp"
#include "StateInputAugmentedLagrangian.hpp"
#include "TestProblemsGeneration.hpp"
#include "iLQR.hpp"

using test_tools::matrix_eigen_conversion::fromEigenVector;
using test_tools::matrix_eigen_conversion::toEigenVector;

class DDPCorrectness : public testing::TestWithParam<unsigned> {
 protected:
  using Scalar = double;
  static constexpr int STATE_DIM = 3;
  static constexpr int INPUT_DIM = 3;
  static constexpr int CONSTRAINT_DIM = 2;
  static constexpr size_t N = 50;
  static constexpr size_t numRandomCases = 3;
  static constexpr int maxFeasibilityRetries = 10;
  static constexpr Scalar timeStep = 1e-2;
  static constexpr Scalar constraintTolerance = 1e-3;
  static constexpr Scalar solutionTolerance = 5e-2;
  static constexpr Scalar inputRelativeTolerance = 5e-2;

  using Transcription_t =
      qp_solver::QpContinuousDynamicsTranscriptionConfig_t<Scalar, STATE_DIM,
                                                           INPUT_DIM, N>;
  using QpProblem_t =
      qp_solver::QpContinuousDynamicsOptimalControlProblem_t<Scalar, STATE_DIM,
                                                             INPUT_DIM, N>;
  using ConstraintConfig_t = ConstraintConfig<
      StateConstraintConfig<ConstraintLayout<>>,
      StateInputConstraintConfig<ConstraintLayout<
          ConstraintGroupLayout<ConstraintTerm<CONSTRAINT_DIM>>,
          ConstraintGroupLayout<>>>,
      FinalStateConstraintConfig<ConstraintLayout<>>>;
  using Descriptor_t =
      iLQRDescriptor<Scalar, Transcription_t, ConstraintConfig_t>;
  using Solver_t = iLQR<Descriptor_t>;
  using AlProblem_t = typename Solver_t::OptimalControlProblem_t;
  using TargetTrajectories_t = typename Solver_t::TargetTrajectories_t;
  using QpTrajectory_t =
      qp_solver::ContinuousTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>;
  using StateVector_t = typename Solver_t::StateVector_t;
  using InputVector_t = typename Solver_t::InputVector_t;
  using Dynamics_t = LinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM>;
  using StateInputCost_t =
      StateInputCost<Scalar, STATE_DIM, INPUT_DIM, static_cast<int>(N + 1)>;
  using FinalCost_t = StateCost<Scalar, STATE_DIM, static_cast<int>(N + 1)>;
  using HardConstraint_t =
      StateInputConstraint<Scalar, STATE_DIM, INPUT_DIM, CONSTRAINT_DIM>;
  using Penalty_t = QuadraticPenalty<Scalar>;
  using AugmentedLagrangian_t =
      StateInputAugmentedLagrangian<Scalar, STATE_DIM, INPUT_DIM,
                                    CONSTRAINT_DIM>;
  using Initializer_t = DefaultInitializer<Scalar, STATE_DIM, INPUT_DIM>;

  struct RandomProblem {
    std::unique_ptr<Dynamics_t> systemPtr;
    std::unique_ptr<StateInputCost_t> qpCostPtr;
    std::unique_ptr<FinalCost_t> qpFinalCostPtr;
    std::unique_ptr<StateInputCost_t> alCostPtr;
    std::unique_ptr<FinalCost_t> alFinalCostPtr;
    std::unique_ptr<HardConstraint_t> hardConstraintPtr;
    std::unique_ptr<Penalty_t> penaltyPtr;
    std::unique_ptr<AugmentedLagrangian_t> augmentedLagrangianPtr;
    std::unique_ptr<QpProblem_t> qpProblemPtr;
    std::unique_ptr<AlProblem_t> alProblemPtr;
    TargetTrajectories_t targetTrajectory;
    QpTrajectory_t nominalTrajectory;
    QpTrajectory_t qpSolution;
    StateVector_t initState;
  };

  DDPCorrectness() {
    std::srand(GetParam());
    randomProblems.reserve(numRandomCases);
    for (size_t caseIndex = 0; caseIndex < numRandomCases; ++caseIndex) {
      randomProblems.emplace_back(createFeasibleRandomProblemWithRetries());
    }
  }

  std::unique_ptr<RandomProblem> createFeasibleRandomProblemWithRetries() {
    for (int retries = 0; retries <= maxFeasibilityRetries; ++retries) {
      auto problem = createFeasibleRandomProblem();
      if (problem) {
        return problem;
      }
      std::cerr << "Random correctness problem was infeasible, retry ...\n";
    }
    throw std::runtime_error("Failed creating feasible correctness problem");
  }

  std::unique_ptr<RandomProblem> createFeasibleRandomProblem() {
    auto problem = std::make_unique<RandomProblem>();

    const auto dynamics =
        test_tools::getRandomDynamics<Scalar, STATE_DIM, INPUT_DIM>();
    problem->systemPtr = test_tools::getiLQRDynamics(dynamics);

    const auto runningCost =
        test_tools::getRandomCost<Scalar, STATE_DIM, INPUT_DIM>();
    const auto finalCost = test_tools::getRandomCost<Scalar, STATE_DIM, 0>();
    problem->qpCostPtr =
        test_tools::getiLQRCost<Scalar, STATE_DIM, INPUT_DIM,
                                static_cast<int>(N + 1)>(runningCost);
    problem->qpFinalCostPtr =
        test_tools::getiLQRStateCost<Scalar, STATE_DIM,
                                     static_cast<int>(N + 1)>(finalCost);
    problem->alCostPtr =
        test_tools::getiLQRCost<Scalar, STATE_DIM, INPUT_DIM,
                                static_cast<int>(N + 1)>(runningCost);
    problem->alFinalCostPtr =
        test_tools::getiLQRStateCost<Scalar, STATE_DIM,
                                     static_cast<int>(N + 1)>(finalCost);

    const auto hardConstraint =
        test_tools::getRandomConstraints<Scalar, STATE_DIM, INPUT_DIM,
                                         CONSTRAINT_DIM>();
    problem->hardConstraintPtr =
        test_tools::getiLQRConstraints<Scalar, STATE_DIM, INPUT_DIM,
                                       CONSTRAINT_DIM>(hardConstraint);

    problem->penaltyPtr = std::make_unique<Penalty_t>(
        typename Penalty_t::Config{Scalar(1e6), Scalar(1.0)});
    problem->augmentedLagrangianPtr = std::make_unique<AugmentedLagrangian_t>(
        problem->hardConstraintPtr.get(), problem->penaltyPtr.get());

    problem->qpProblemPtr = std::make_unique<QpProblem_t>();
    problem->qpProblemPtr->dynamicsPtr = problem->systemPtr.get();
    problem->qpProblemPtr->cost.add(*problem->qpCostPtr);
    problem->qpProblemPtr->finalCost.add(*problem->qpFinalCostPtr);

    problem->alProblemPtr = std::make_unique<AlProblem_t>();
    problem->alProblemPtr->dynamicsPtr = problem->systemPtr.get();
    problem->alProblemPtr->cost.add(*problem->alCostPtr);
    problem->alProblemPtr->finalCost.add(*problem->alFinalCostPtr);
    problem->alProblemPtr->equalityLagrangian.template set<0>(
        problem->augmentedLagrangianPtr.get());

    setRandomTargetTrajectory(*problem);
    setDynamicallyConsistentNominalTrajectory(*problem);
    problem->initState = StateVector_t::Random();

    const auto lqApproximation = qp_solver::getLinearQuadraticApproximation(
        *problem->qpProblemPtr, problem->targetTrajectory,
        problem->nominalTrajectory, *problem->hardConstraintPtr);
    const Eigen::Vector<Scalar, STATE_DIM> dx0 =
        toEigenVector(problem->initState) -
        problem->nominalTrajectory.stateTrajectory.front();
    const auto denseQp = qp_solver::getDenseQp(lqApproximation, dx0);

    if (!test_tools::isQpFeasible(denseQp.first, denseQp.second)) {
      return nullptr;
    }

    problem->nominalTrajectory =
        getFeasibleTrajectory(denseQp.second, problem->nominalTrajectory);

    try {
      problem->qpSolution = solveHardConstrainedQp(*problem);
    } catch (const std::runtime_error&) {
      return nullptr;
    }
    return problem;
  }

  /** 修改给定轨迹以满足约束。 */
  QpTrajectory_t getFeasibleTrajectory(
      const qp_solver::QpDenseConstraintApproximation<Scalar>& qpConstraints,
      const QpTrajectory_t& trajectory) const {
    const auto& A = qpConstraints.dfdx;  // A w + b = 0，A 必须行满秩，使得
                                         // (A A') 可逆。
    const auto& b =
        qpConstraints.f;  // b = [x0; e[0]; b[0]; ... e[N-1]; b[N-1]; e[N]]

    /* 通过求解以下问题寻找满足约束的轨迹修正 w：
     *   min  1/2 w' w。
     *   s.t. A w + b = 0  */
    const Eigen::Vector<Scalar, Eigen::Dynamic> w =
        -A.transpose() * (A * A.transpose()).inverse() * b;

    // 使轨迹可行。
    auto feasibleTrajectory = trajectory;
    Eigen::Index nextIndex = 0;
    const auto numStages = feasibleTrajectory.inputTrajectory.size();
    for (size_t k = 0; k < numStages; k++) {
      const auto nx = feasibleTrajectory.stateTrajectory[k].size();
      const auto nu = feasibleTrajectory.inputTrajectory[k].size();
      feasibleTrajectory.stateTrajectory[k] +=
          w.segment(nextIndex, nx);  // dx[k]
      feasibleTrajectory.inputTrajectory[k] +=
          w.segment(nextIndex + nx, nu);  // du[k]
      nextIndex += nx + nu;
    }
    feasibleTrajectory.stateTrajectory[numStages] += w.segment(
        nextIndex,
        feasibleTrajectory.stateTrajectory[numStages].size());  // dx[N]

    return feasibleTrajectory;
  }

  void setRandomTargetTrajectory(RandomProblem& problem) {
    for (size_t k = 0; k < N + 1; ++k) {
      problem.targetTrajectory.timeTrajectory[k] =
          static_cast<Scalar>(k) * timeStep;
      problem.targetTrajectory.stateTrajectory[k] = StateVector_t::Random();
      // The target input array has N+1 slots to match the solver trajectory
      // convention, even though the QP trajectory uses only N inputs.
      problem.targetTrajectory.inputTrajectory[k] = InputVector_t::Random();
    }
  }

  void setDynamicallyConsistentNominalTrajectory(RandomProblem& problem) {
    EK2DynamicsDiscretizer<Scalar, STATE_DIM, INPUT_DIM> discretizer;
    problem.nominalTrajectory.stateTrajectory[0] =
        toEigenVector(StateVector_t::Random());
    for (size_t k = 0; k < N; ++k) {
      problem.nominalTrajectory.timeTrajectory[k] =
          static_cast<Scalar>(k) * timeStep;
      problem.nominalTrajectory.inputTrajectory[k] =
          toEigenVector(InputVector_t::Random());
      const auto nextState = discretizer.discretize(
          *problem.systemPtr, problem.nominalTrajectory.timeTrajectory[k],
          fromEigenVector(problem.nominalTrajectory.stateTrajectory[k]),
          fromEigenVector(problem.nominalTrajectory.inputTrajectory[k]),
          timeStep);
      problem.nominalTrajectory.stateTrajectory[k + 1] =
          toEigenVector(nextState);
    }
    problem.nominalTrajectory.timeTrajectory[N] =
        static_cast<Scalar>(N) * timeStep;
  }

  QpTrajectory_t solveHardConstrainedQp(const RandomProblem& problem) const {
    const auto lqApproximation = qp_solver::getLinearQuadraticApproximation(
        *problem.qpProblemPtr, problem.targetTrajectory,
        problem.nominalTrajectory, *problem.hardConstraintPtr);
    const Eigen::Vector<Scalar, STATE_DIM> dx0 =
        toEigenVector(problem.initState) -
        problem.nominalTrajectory.stateTrajectory.front();
    const auto deltaTrajectories =
        qp_solver::solveLinearQuadraticProblem(lqApproximation, dx0);

    QpTrajectory_t deltaSolution;
    deltaSolution.timeTrajectory = problem.nominalTrajectory.timeTrajectory;
    for (size_t k = 0; k < N + 1; ++k) {
      deltaSolution.stateTrajectory[k] = deltaTrajectories.first[k];
    }
    for (size_t k = 0; k < N; ++k) {
      deltaSolution.inputTrajectory[k] = deltaTrajectories.second[k];
    }
    return problem.nominalTrajectory + deltaSolution;
  }

  DDPSettings<Scalar> getSettings() const {
    DDPSettings<Scalar> ddpSettings;
    ddpSettings.timeStep_ = timeStep;
    ddpSettings.maxNumIterations_ = 30;
    ddpSettings.minRelCost_ = 1e-9;
    ddpSettings.constraintTolerance_ = 1e-6;
    ddpSettings.strategy_ = SearchStrategyType::LINE_SEARCH;
    ddpSettings.lineSearch_.minStepLength = 1e-2;
    ddpSettings.lineSearch_.maxStepLength = 1.0;
    ddpSettings.lineSearch_.hessianCorrectionMultiple = 1e-6;
    return ddpSettings;
  }

  Scalar maxConstraintViolation(
      const RandomProblem& problem,
      const typename Solver_t::PrimalSolution_t& solution) const {
    Scalar maxViolation = Scalar(0.0);
    for (size_t k = 0; k < N; ++k) {
      maxViolation =
          std::max(maxViolation, problem.hardConstraintPtr
                                     ->getValue(solution.timeTrajectory_[k],
                                                solution.stateTrajectory_[k],
                                                solution.inputTrajectory_[k])
                                     .template lpNorm<matrix::Infinity>());
    }
    return maxViolation;
  }

  Scalar maxConstraintViolation(const RandomProblem& problem,
                                const QpTrajectory_t& solution) const {
    Scalar maxViolation = Scalar(0.0);
    for (size_t k = 0; k < N; ++k) {
      maxViolation =
          std::max(maxViolation,
                   problem.hardConstraintPtr
                       ->getValue(solution.timeTrajectory[k],
                                  fromEigenVector(solution.stateTrajectory[k]),
                                  fromEigenVector(solution.inputTrajectory[k]))
                       .template lpNorm<matrix::Infinity>());
    }
    return maxViolation;
  }

  std::vector<std::unique_ptr<RandomProblem>> randomProblems;
  Initializer_t initializer;
};

TEST_P(DDPCorrectness, HardQpSolutionSatisfiesLinearConstraints) {
  for (size_t caseIndex = 0; caseIndex < numRandomCases; ++caseIndex) {
    const auto& problem = *randomProblems[caseIndex];
    EXPECT_LT(maxConstraintViolation(problem, problem.qpSolution),
              constraintTolerance)
        << "seed: " << GetParam() << ", case index: " << caseIndex;
  }
}

TEST_P(DDPCorrectness, ILQRWithAugmentedLagrangianMatchesHardQpSolution) {
  for (size_t caseIndex = 0; caseIndex < numRandomCases; ++caseIndex) {
    const auto& problem = *randomProblems[caseIndex];

    Solver_t solver(getSettings(), *problem.alProblemPtr, &initializer);
    solver.setDesireTrajectory(problem.targetTrajectory.timeTrajectory,
                               problem.targetTrajectory.stateTrajectory,
                               problem.targetTrajectory.inputTrajectory);

    ASSERT_NO_THROW(solver.run(Scalar(0.0), problem.initState))
        << "seed: " << GetParam() << ", case index: " << caseIndex;
    const auto& solution = solver.primalSolution();

    EXPECT_LT(maxConstraintViolation(problem, solution), constraintTolerance)
        << "seed: " << GetParam() << ", case index: " << caseIndex;
    EXPECT_LT((toEigenVector(solution.stateTrajectory_.back()) -
               problem.qpSolution.stateTrajectory.back())
                  .norm(),
              solutionTolerance)
        << "seed: " << GetParam() << ", case index: " << caseIndex
        << ", Hard QP final state: "
        << problem.qpSolution.stateTrajectory.back().transpose()
        << ", iLQR final state: "
        << toEigenVector(solution.stateTrajectory_.back()).transpose();
    const Scalar initialInputRelativeError =
        (toEigenVector(solution.inputTrajectory_.front()) -
         problem.qpSolution.inputTrajectory.front())
            .norm() /
        std::max(problem.qpSolution.inputTrajectory.front().norm(),
                 Scalar(1.0));
    EXPECT_LT(initialInputRelativeError, inputRelativeTolerance)
        << "seed: " << GetParam() << ", case index: " << caseIndex
        << ", Hard QP initial input: "
        << problem.qpSolution.inputTrajectory.front().transpose()
        << ", iLQR initial input: "
        << toEigenVector(solution.inputTrajectory_.front()).transpose();
  }
}

INSTANTIATE_TEST_SUITE_P(
    RandomSeeds, DDPCorrectness, testing::Values(0U, 1U, 2U),
    [](const testing::TestParamInfo<DDPCorrectness::ParamType>& info) {
      return "Seed" + std::to_string(info.param);
    });
