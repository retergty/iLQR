#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

#include "DefaultInitializer.hpp"
#include "DiscreteLinearSystemDynamics.hpp"
#include "QpDiscreteTranscription.hpp"
#include "QpSolver.hpp"
#include "QuadraticPenalty.hpp"
#include "StateInputAugmentedLagrangian.hpp"
#include "TestProblemsGeneration.hpp"
#include "iLQR.hpp"


class DiscreteDDPCorrectness : public testing::TestWithParam<unsigned> {
 protected:
  using Scalar = double;
  static constexpr int STATE_DIM = 3;
  static constexpr int INPUT_DIM = 3;
  static constexpr int CONSTRAINT_DIM = 2;
  static constexpr size_t N = 20;
  static constexpr size_t numRandomCases = 2;
  static constexpr int maxFeasibilityRetries = 10;
  static constexpr Scalar timeStep = 1e-2;
  static constexpr Scalar constraintTolerance = 1e-3;
  static constexpr Scalar solutionTolerance = 5e-2;
  static constexpr Scalar inputRelativeTolerance = 5e-2;

  using Transcription_t =
      qp_solver::QpDiscreteDynamicsTranscriptionConfig_t<Scalar, STATE_DIM,
                                                         INPUT_DIM, N>;
  using QpProblem_t =
      qp_solver::QpDiscreteDynamicsOptimalControlProblem_t<Scalar, STATE_DIM,
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
  using StateVector_t = Vector<Scalar, STATE_DIM>;
  using InputVector_t = Vector<Scalar, INPUT_DIM>;
  using Dynamics_t = DiscreteLinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM>;
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

  DiscreteDDPCorrectness() {
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
      std::cerr << "Discrete random correctness problem was infeasible, retry "
                   "...\n";
    }
    throw std::runtime_error(
        "Failed creating feasible discrete correctness problem");
  }

  std::unique_ptr<RandomProblem> createFeasibleRandomProblem() {
    auto problem = std::make_unique<RandomProblem>();

    Matrix<Scalar, STATE_DIM, STATE_DIM> A =
        Matrix<Scalar, STATE_DIM, STATE_DIM>::Identity();
    A += Scalar(0.05) * Matrix<Scalar, STATE_DIM, STATE_DIM>::Random();
    Matrix<Scalar, STATE_DIM, INPUT_DIM> B =
        Scalar(0.1) * Matrix<Scalar, STATE_DIM, INPUT_DIM>::Random();
    problem->systemPtr = std::make_unique<Dynamics_t>(A, B);

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
    const StateVector_t dx0 =
        problem->initState - problem->nominalTrajectory.stateTrajectory.front();
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

  QpTrajectory_t getFeasibleTrajectory(
      const VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic,
                                              Eigen::Dynamic, 0>& qpConstraints,
      const QpTrajectory_t& trajectory) const {
    const auto& A = qpConstraints.dfdx;
    const auto& b = qpConstraints.f;
    const Vector<Scalar, Eigen::Dynamic> w =
        -A.transpose() * (A * A.transpose()).inverse() * b;

    auto feasibleTrajectory = trajectory;
    Eigen::Index nextIndex = 0;
    const auto numStages = feasibleTrajectory.inputTrajectory.size();
    for (size_t k = 0; k < numStages; k++) {
      const auto nx = feasibleTrajectory.stateTrajectory[k].size();
      const auto nu = feasibleTrajectory.inputTrajectory[k].size();
      feasibleTrajectory.stateTrajectory[k] += w.segment(nextIndex, nx);
      feasibleTrajectory.inputTrajectory[k] += w.segment(nextIndex + nx, nu);
      nextIndex += nx + nu;
    }
    feasibleTrajectory.stateTrajectory[numStages] += w.segment(
        nextIndex, feasibleTrajectory.stateTrajectory[numStages].size());
    return feasibleTrajectory;
  }

  void setRandomTargetTrajectory(RandomProblem& problem) {
    for (size_t k = 0; k < N + 1; ++k) {
      problem.targetTrajectory.timeTrajectory[k] =
          static_cast<Scalar>(k) * timeStep;
      problem.targetTrajectory.stateTrajectory[k] = StateVector_t::Random();
      problem.targetTrajectory.inputTrajectory[k] = InputVector_t::Random();
    }
  }

  void setDynamicallyConsistentNominalTrajectory(RandomProblem& problem) {
    problem.nominalTrajectory.stateTrajectory[0] = StateVector_t::Random();
    for (size_t k = 0; k < N; ++k) {
      problem.nominalTrajectory.timeTrajectory[k] =
          static_cast<Scalar>(k) * timeStep;
      problem.nominalTrajectory.inputTrajectory[k] = InputVector_t::Random();
      problem.nominalTrajectory.stateTrajectory[k + 1] =
          problem.systemPtr->computeMap(
              problem.nominalTrajectory.timeTrajectory[k],
              problem.nominalTrajectory.stateTrajectory[k],
              problem.nominalTrajectory.inputTrajectory[k], timeStep);
    }
    problem.nominalTrajectory.timeTrajectory[N] =
        static_cast<Scalar>(N) * timeStep;
  }

  QpTrajectory_t solveHardConstrainedQp(const RandomProblem& problem) const {
    const auto lqApproximation = qp_solver::getLinearQuadraticApproximation(
        *problem.qpProblemPtr, problem.targetTrajectory,
        problem.nominalTrajectory, *problem.hardConstraintPtr);
    const StateVector_t dx0 =
        problem.initState - problem.nominalTrajectory.stateTrajectory.front();
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
                                     .template lpNorm<Eigen::Infinity>());
    }
    return maxViolation;
  }

  Scalar maxConstraintViolation(const RandomProblem& problem,
                                const QpTrajectory_t& solution) const {
    Scalar maxViolation = Scalar(0.0);
    for (size_t k = 0; k < N; ++k) {
      maxViolation =
          std::max(maxViolation, problem.hardConstraintPtr
                                     ->getValue(solution.timeTrajectory[k],
                                                solution.stateTrajectory[k],
                                                solution.inputTrajectory[k])
                                     .template lpNorm<Eigen::Infinity>());
    }
    return maxViolation;
  }

  std::vector<std::unique_ptr<RandomProblem>> randomProblems;
  Initializer_t initializer;
};

TEST_P(DiscreteDDPCorrectness, HardQpSolutionSatisfiesLinearConstraints) {
  for (size_t caseIndex = 0; caseIndex < numRandomCases; ++caseIndex) {
    const auto& problem = *randomProblems[caseIndex];
    EXPECT_LT(maxConstraintViolation(problem, problem.qpSolution),
              constraintTolerance)
        << "seed: " << GetParam() << ", case index: " << caseIndex;
  }
}

TEST_P(DiscreteDDPCorrectness,
       ILQRWithAugmentedLagrangianMatchesHardQpSolution) {
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
    EXPECT_LT((solution.stateTrajectory_.back() -
               problem.qpSolution.stateTrajectory.back())
                  .norm(),
              solutionTolerance)
        << "seed: " << GetParam() << ", case index: " << caseIndex;
    const Scalar initialInputRelativeError =
        (solution.inputTrajectory_.front() -
         problem.qpSolution.inputTrajectory.front())
            .norm() /
        std::max(problem.qpSolution.inputTrajectory.front().norm(),
                 Scalar(1.0));
    EXPECT_LT(initialInputRelativeError, inputRelativeTolerance)
        << "seed: " << GetParam() << ", case index: " << caseIndex;
  }
}

INSTANTIATE_TEST_SUITE_P(
    RandomSeeds, DiscreteDDPCorrectness, testing::Values(0U, 1U),
    [](const testing::TestParamInfo<DiscreteDDPCorrectness::ParamType>& info) {
      return "Seed" + std::to_string(info.param);
    });
