/******************************************************************************
Copyright (c) 2017, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <gtest/gtest.h>

#include <cstdlib>
#include <ctime>
#include <memory>
#include <string>

#include "DefaultInitializer.hpp"
#include "EXP0.hpp"
#include "Ocs2QpSolver.hpp"
#include "iLQR.hpp"

class Exp0 : public testing::Test {
 protected:
  using Scalar = double;
  static constexpr int STATE_DIM = exp0::STATE_DIM;
  static constexpr int INPUT_DIM = exp0::INPUT_DIM;
  static constexpr Scalar timeStep = 1e-2;
  static constexpr Scalar minRelCost = 1e-3;
  static constexpr size_t PredictLength = 200;

  using Descriptor_t =
      iLQRDescriptor<Scalar,
                     TranscriptionConfig<Dimensions<STATE_DIM, INPUT_DIM>,
                                         Horizon<PredictLength>>>;
  using Solver_t = iLQR<Descriptor_t>;
  using Problem_t = typename Solver_t::OptimalControlProblem_t;
  using PrimalSolution_t = typename Solver_t::PrimalSolution_t;
  using DualSolution_t = typename Solver_t::DualSolution_t;
  using ProblemMetrics_t = typename Solver_t::ProblemMetrics_t;
  using QpTrajectory_t =
      qp_solver::ContinuousTrajectory<Scalar, STATE_DIM, INPUT_DIM,
                                      PredictLength>;
  using Initializer_t = DefaultInitializer<Scalar, STATE_DIM, INPUT_DIM>;
  using StateVector_t = Vector<Scalar, STATE_DIM>;
  using InputVector_t = Vector<Scalar, INPUT_DIM>;
  using RolloutSettings_t = RolloutSettings<Scalar>;
  using DDPSettings_t = DDPSettings<Scalar>;

  Exp0()
      : problem(exp0::createExp0Problem<Scalar, PredictLength>()),
        initializerPtr(std::make_unique<Initializer_t>()) {
    setConstantReferenceTrajectory();
  }

  // rollout settings
  RolloutSettings_t rolloutSettings() const {
    RolloutSettings_t rolloutSettings;
    rolloutSettings.timeStep = timeStep;
    return rolloutSettings;
  };

  DDPSettings_t getSettings(SearchStrategyType strategy,
                            bool useFeedbackPolicy = true) const {
    DDPSettings_t ddpSettings;
    ddpSettings.timeStep_ = timeStep;
    ddpSettings.maxNumIterations_ = 30;
    ddpSettings.minRelCost_ = minRelCost;
    ddpSettings.useFeedbackPolicy_ = useFeedbackPolicy;
    ddpSettings.strategy_ = strategy;
    ddpSettings.lineSearch_.minStepLength = 0.0001;
    return ddpSettings;
  }

  std::string getTestName(const DDPSettings_t& ddpSettings) const {
    std::string testName;
    testName += "EXP0 Test { ";
    testName += "Algorithm: iLQR,  ";
    testName +=
        "Strategy: " + searchStrategyToString(ddpSettings.strategy_) + " }";
    return testName;
  }

  void setConstantReferenceTrajectory() {
    const StateVector_t targetState = exp0::getExp0TargetState<Scalar>();
    const InputVector_t targetInput = exp0::getExp0TargetInput<Scalar>();

    for (size_t k = 0; k < PredictLength + 1; ++k) {
      problem.timeTrajectory[k] = static_cast<Scalar>(k) * timeStep;
      problem.stateTrajectory[k] = targetState;
      problem.inputTrajectory[k] = targetInput;
    }
  }

  QpTrajectory_t toQpTrajectory(const PrimalSolution_t& solution) const {
    QpTrajectory_t trajectory;
    for (size_t k = 0; k < PredictLength + 1; ++k) {
      trajectory.timeTrajectory[k] = solution.timeTrajectory_[k];
      trajectory.stateTrajectory[k] = solution.stateTrajectory_[k];
    }
    for (size_t k = 0; k < PredictLength; ++k) {
      trajectory.inputTrajectory[k] = solution.inputTrajectory_[k];
    }
    return trajectory;
  }

  PrimalSolution_t toPrimalSolution(const QpTrajectory_t& trajectory) const {
    PrimalSolution_t solution;
    for (size_t k = 0; k < PredictLength + 1; ++k) {
      solution.timeTrajectory_[k] = trajectory.timeTrajectory[k];
      solution.stateTrajectory_[k] = trajectory.stateTrajectory[k];
    }
    for (size_t k = 0; k < PredictLength; ++k) {
      solution.inputTrajectory_[k] = trajectory.inputTrajectory[k];
    }
    solution.inputTrajectory_[PredictLength] =
        trajectory.inputTrajectory[PredictLength - 1];
    return solution;
  }

  typename Solver_t::PerformanceIndex_t getPerformanceIndex(
      Problem_t& targetProblem, const PrimalSolution_t& solution) const {
    DualSolution_t dualSolution;
    ProblemMetrics_t problemMetrics;
    Solver_t::computeRolloutMetrics(targetProblem, solution, dualSolution,
                                    problemMetrics);
    return Solver_t::computeRolloutPerformanceIndex(solution.timeTrajectory_,
                                                    problemMetrics);
  }

  static std::string searchStrategyToString(SearchStrategyType strategy) {
    switch (strategy) {
      case SearchStrategyType::LINE_SEARCH:
        return "LINE_SEARCH";
      case SearchStrategyType::LEVENBERG_MARQUARDT:
        return "LEVENBERG_MARQUARDT";
      default:
        return "UNKNOWN";
    }
  }

  const Scalar startTime = 0.0;
  const Scalar finalTime = 2.0;
  const StateVector_t initState = (StateVector_t() << 0.0, 2.0).finished();

  Problem_t& problem;
  std::unique_ptr<Initializer_t> initializerPtr;
};

constexpr int Exp0::STATE_DIM;
constexpr int Exp0::INPUT_DIM;
constexpr Exp0::Scalar Exp0::timeStep;
constexpr Exp0::Scalar Exp0::minRelCost;
constexpr size_t Exp0::PredictLength;

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, ddp_feedback_policy) {
  // ddp settings
  auto ddpSettings = getSettings(SearchStrategyType::LINE_SEARCH, true);

  // dynamics and rollout
  exp0::EXP0_Sys1<Scalar> systemDynamics;

  // instantiate
  Solver_t ddp(ddpSettings, &systemDynamics, initializerPtr.get());
  exp0::EXP0_Cost<Scalar, static_cast<int>(PredictLength + 1)> cost;
  exp0::EXP0_FinalCost<Scalar, static_cast<int>(PredictLength + 1)> finalCost;
  ddp.optimalControlProblem_.cost.add(cost);
  ddp.optimalControlProblem_.finalCost.add(finalCost);
  ddp.setDesireTrajectory(problem.timeTrajectory, problem.stateTrajectory,
                          problem.inputTrajectory);

  // run ddp
  EXPECT_NO_THROW(ddp.run(startTime, initState));

  // get solution
  const auto& solution = ddp.primalSolution();
  const auto& ctrl = solution.controller_;

  EXPECT_EQ(ctrl.getType(), ControllerType::LINEAR)
      << "MESSAGE: iLQR solution does not contain a linear feedback policy!";
  EXPECT_DOUBLE_EQ(ctrl.timeStamp_.back(), finalTime)
      << "MESSAGE: iLQR failed in policy final time of controller!";
  EXPECT_DOUBLE_EQ(solution.timeTrajectory_.back(), finalTime)
      << "MESSAGE: iLQR failed in policy final time of trajectory!";
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, ddp_feedforward_policy) {
  // ddp settings
  auto ddpSettings = getSettings(SearchStrategyType::LINE_SEARCH, false);

  // dynamics and rollout
  exp0::EXP0_Sys1<Scalar> systemDynamics;

  // instantiate
  Solver_t ddp(ddpSettings, &systemDynamics, initializerPtr.get());
  exp0::EXP0_Cost<Scalar, static_cast<int>(PredictLength + 1)> cost;
  exp0::EXP0_FinalCost<Scalar, static_cast<int>(PredictLength + 1)> finalCost;
  ddp.optimalControlProblem_.cost.add(cost);
  ddp.optimalControlProblem_.finalCost.add(finalCost);
  ddp.setDesireTrajectory(problem.timeTrajectory, problem.stateTrajectory,
                          problem.inputTrajectory);

  // run ddp
  EXPECT_NO_THROW(ddp.run(startTime, initState));

  // get solution
  const auto& solution = ddp.primalSolution();
  const auto& ctrl = solution.controller_;

  EXPECT_EQ(ctrl.getType(), ControllerType::LINEAR)
      << "MESSAGE: iLQR solution does not contain a controller with "
         "feedforward terms!";
  EXPECT_TRUE(ctrl.biasArray_.back().allFinite())
      << "MESSAGE: iLQR feedforward policy contains invalid values!";
  EXPECT_DOUBLE_EQ(ctrl.timeStamp_.back(), finalTime)
      << "MESSAGE: iLQR failed in policy final time of controller!";
  EXPECT_DOUBLE_EQ(solution.timeTrajectory_.back(), finalTime)
      << "MESSAGE: iLQR failed in policy final time of trajectory!";
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, ddp_moving_horizon) {
  // ddp settings
  const auto ddpSettings = getSettings(SearchStrategyType::LINE_SEARCH, true);

  // dynamics and rollout
  exp0::EXP0_Sys1<Scalar> systemDynamics;

  // instantiate
  Solver_t ddp(ddpSettings, &systemDynamics, initializerPtr.get());
  exp0::EXP0_Cost<Scalar, static_cast<int>(PredictLength + 1)> cost;
  exp0::EXP0_FinalCost<Scalar, static_cast<int>(PredictLength + 1)> finalCost;
  ddp.optimalControlProblem_.cost.add(cost);
  ddp.optimalControlProblem_.finalCost.add(finalCost);
  ddp.setDesireTrajectory(problem.timeTrajectory, problem.stateTrajectory,
                          problem.inputTrajectory);

  const auto expectSolutionEndsAt = [&ddp](Scalar expectedFinalTime) {
    const auto& solution = ddp.primalSolution();
    EXPECT_DOUBLE_EQ(solution.controller_.timeStamp_.back(), expectedFinalTime)
        << "MESSAGE: iLQR failed in policy final time of controller!";
    EXPECT_DOUBLE_EQ(solution.timeTrajectory_.back(), expectedFinalTime)
        << "MESSAGE: iLQR failed in policy final time of trajectory!";
  };

  // run iLQR (no active event)
  Scalar movingStartTime = 0.2;
  EXPECT_NO_THROW(ddp.run(movingStartTime, initState));
  expectSolutionEndsAt(movingStartTime +
                       static_cast<Scalar>(PredictLength) * timeStep);

  // move the time horizon forward with overlap
  movingStartTime = 0.6;
  EXPECT_NO_THROW(ddp.run(movingStartTime, initState));
  expectSolutionEndsAt(movingStartTime +
                       static_cast<Scalar>(PredictLength) * timeStep);

  // move the time horizon forward with partial overlap
  movingStartTime = 1.1;
  EXPECT_NO_THROW(ddp.run(movingStartTime, initState));
  expectSolutionEndsAt(movingStartTime +
                       static_cast<Scalar>(PredictLength) * timeStep);

  // move the time horizon (no overlap)
  movingStartTime = 1.6;
  EXPECT_NO_THROW(ddp.run(movingStartTime, initState));
  expectSolutionEndsAt(movingStartTime +
                       static_cast<Scalar>(PredictLength) * timeStep);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, qp_solver_matches_ilqr_solution) {
  auto ddpSettings = getSettings(SearchStrategyType::LINE_SEARCH, true);
  ddpSettings.maxNumIterations_ = 50;
  ddpSettings.minRelCost_ = 1e-9;

  exp0::EXP0_Sys1<Scalar> systemDynamics;

  Solver_t ddp(ddpSettings, &systemDynamics, initializerPtr.get());
  exp0::EXP0_Cost<Scalar, static_cast<int>(PredictLength + 1)> cost;
  exp0::EXP0_FinalCost<Scalar, static_cast<int>(PredictLength + 1)> finalCost;
  ddp.optimalControlProblem_.cost.add(cost);
  ddp.optimalControlProblem_.finalCost.add(finalCost);
  ddp.setDesireTrajectory(problem.timeTrajectory, problem.stateTrajectory,
                          problem.inputTrajectory);

  ASSERT_NO_THROW(ddp.run(startTime, initState));
  const auto& ilqrSolution = ddp.primalSolution();
  const auto ilqrPerformance = ddp.performanceIndex();

  const auto qpSolution = qp_solver::solveLinearQuadraticOptimalControlProblem(
      ddp.optimalControlProblem_, toQpTrajectory(ilqrSolution), initState);
  const auto qpPrimalSolution = toPrimalSolution(qpSolution);
  const auto qpPerformance =
      getPerformanceIndex(ddp.optimalControlProblem_, qpPrimalSolution);

  EXPECT_NEAR(qpPerformance.cost, ilqrPerformance.cost, 1e-1)
      << "QP cost: " << qpPerformance.cost
      << ", iLQR cost: " << ilqrPerformance.cost;
  EXPECT_LT(
      (qpSolution.stateTrajectory.back() - ilqrSolution.stateTrajectory_.back())
          .norm(),
      1e-1)
      << "QP final state: " << qpSolution.stateTrajectory.back().transpose()
      << ", iLQR final state: "
      << ilqrSolution.stateTrajectory_.back().transpose();
  EXPECT_LT((qpSolution.inputTrajectory.front() -
             ilqrSolution.inputTrajectory_.front())
                .norm(),
            1e-1)
      << "QP initial input: " << qpSolution.inputTrajectory.front().transpose()
      << ", iLQR initial input: "
      << ilqrSolution.inputTrajectory_.front().transpose();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, ddp_q_function) {
  // ddp settings
  auto ddpSettings = getSettings(SearchStrategyType::LINE_SEARCH, true);
  ddpSettings.maxNumIterations_ = 50;
  ddpSettings.minRelCost_ = 1e-9;  // to allow more iterations that the effect
                                   // of final linesearch is negligible

  // dynamics and rollout
  exp0::EXP0_Sys1<Scalar> systemDynamics;

  // instantiate
  Solver_t ddp(ddpSettings, &systemDynamics, initializerPtr.get());
  exp0::EXP0_Cost<Scalar, static_cast<int>(PredictLength + 1)> cost;
  exp0::EXP0_FinalCost<Scalar, static_cast<int>(PredictLength + 1)> finalCost;
  ddp.optimalControlProblem_.cost.add(cost);
  ddp.optimalControlProblem_.finalCost.add(finalCost);
  ddp.setDesireTrajectory(problem.timeTrajectory, problem.stateTrajectory,
                          problem.inputTrajectory);

  // run ddp
  ddp.run(startTime, initState);
  // get solution
  const auto& solution = ddp.primalSolution();
  const auto& controller = solution.controller_;

  // define precision for tests
  constexpr Scalar precision = 1e-3;

  const size_t timeIndex = 0;

  // get Q function at current solution
  // expected outcome: true, because the current solution should be optimal
  const Scalar time = solution.timeTrajectory_[timeIndex];
  StateVector_t state = solution.stateTrajectory_[timeIndex];
  InputVector_t input = controller.computeInput(time, state);
  auto qFunction = ddp.getQFunction(timeIndex, state, input);
  const InputVector_t dQdu1a = qFunction.dfdu;
  EXPECT_TRUE(dQdu1a.isZero(precision))
      << "MESSAGE for test 1a: Derivative of Q function w.r.t. to u is not "
         "zero: "
      << dQdu1a.transpose();

  // evaluate Q function at different state (but using feedback policy)
  // expected outcome: true, because for a linear system the LQ approximation of
  // Q is exact and the linear feedback policy is globally optimal
  StateVector_t queryState = StateVector_t::Random();
  InputVector_t queryInput = controller.computeInput(time, queryState);
  const InputVector_t dQdu1b = qFunction.dfdux * (queryState - state) +
                               qFunction.dfduu * (queryInput - input) +
                               qFunction.dfdu;
  EXPECT_TRUE(dQdu1b.isZero(precision))
      << "MESSAGE for test 1b: Derivative of Q function w.r.t. to u is not "
         "zero: "
      << dQdu1b.transpose();

  // evaluate Q function at different input
  // expected outcome: false, because for a linear system the LQ approximation
  // of Q is exact and a random input is not optimal
  queryState = solution.stateTrajectory_[timeIndex];
  queryInput = InputVector_t::Random();
  const InputVector_t dQdu1c = qFunction.dfdux * (queryState - state) +
                               qFunction.dfduu * (queryInput - input) +
                               qFunction.dfdu;
  EXPECT_FALSE(dQdu1c.isZero(precision))
      << "MESSAGE for test 1c: Derivative of Q function w.r.t. to u is zero: "
      << dQdu1c.transpose();

  // get Q function at different state (but using feedback policy)
  // expected outcome: true, because for a linear system the linear feedback
  // policy is globally optimal
  state = StateVector_t::Random();
  input = controller.computeInput(time, state);
  qFunction = ddp.getQFunction(timeIndex, state, input);
  const InputVector_t dQdu2 = qFunction.dfdu;
  EXPECT_TRUE(dQdu2.isZero(precision)) << "MESSAGE for test 2: Derivative of Q "
                                          "function w.r.t. to u is not zero: "
                                       << dQdu2.transpose();

  // get Q function at different input
  // expected outcome: false, because a random input is not optimal
  state = solution.stateTrajectory_[timeIndex];
  input = InputVector_t::Random();
  qFunction = ddp.getQFunction(timeIndex, state, input);
  const InputVector_t dQdu3 = qFunction.dfdu;
  EXPECT_FALSE(dQdu3.isZero(precision))
      << "MESSAGE for test 3: Derivative of Q function w.r.t. to u is zero: "
      << dQdu3.transpose();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/* Add parameterized test suite */
class Exp0Param : public Exp0,
                  public testing::WithParamInterface<SearchStrategyType> {
 protected:
  SearchStrategyType getSearchStrategy() const { return GetParam(); }
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_P(Exp0Param, ILQR) {
  // ddp settings
  const auto ddpSettings = getSettings(getSearchStrategy(), true);

  // dynamics and rollout
  exp0::EXP0_Sys1<Scalar> systemDynamics;

  // instantiate
  Solver_t ddp(ddpSettings, &systemDynamics, initializerPtr.get());
  exp0::EXP0_Cost<Scalar, static_cast<int>(PredictLength + 1)> cost;
  exp0::EXP0_FinalCost<Scalar, static_cast<int>(PredictLength + 1)> finalCost;
  ddp.optimalControlProblem_.cost.add(cost);
  ddp.optimalControlProblem_.finalCost.add(finalCost);
  ddp.setDesireTrajectory(problem.timeTrajectory, problem.stateTrajectory,
                          problem.inputTrajectory);

  // run ddp
  EXPECT_NO_THROW(ddp.run(startTime, initState));

  // get solution
  const auto& solution = ddp.primalSolution();
  const auto& ctrl = solution.controller_;

  EXPECT_EQ(ctrl.getType(), ControllerType::LINEAR)
      << "MESSAGE: " << getTestName(ddpSettings)
      << ": iLQR solution does not contain a linear feedback policy!";
  EXPECT_DOUBLE_EQ(ctrl.timeStamp_.back(), finalTime)
      << "MESSAGE: " << getTestName(ddpSettings)
      << ": failed in policy final time of controller!";
  EXPECT_DOUBLE_EQ(solution.timeTrajectory_.back(), finalTime)
      << "MESSAGE: " << getTestName(ddpSettings)
      << ": failed in policy final time of trajectory!";
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
INSTANTIATE_TEST_SUITE_P(
    Exp0ParamCase, Exp0Param,
    testing::ValuesIn({SearchStrategyType::LINE_SEARCH,
                       SearchStrategyType::LEVENBERG_MARQUARDT}),
    [](const testing::TestParamInfo<Exp0Param::ParamType>& info) {
      /* returns test name for gtest summary */
      switch (info.param) {
        case SearchStrategyType::LINE_SEARCH:
          return std::string("LINE_SEARCH");
        case SearchStrategyType::LEVENBERG_MARQUARDT:
          return std::string("LEVENBERG_MARQUARDT");
        default:
          return std::string("UNKNOWN");
      }
    });
