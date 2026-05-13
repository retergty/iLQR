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
#include <iostream>
#include <memory>
#include <string>

#include "DefaultInitializer.hpp"
#include "EXP0.hpp"
#include "iLQR.hpp"

class Exp0 : public testing::Test {
protected:
  using Scalar = double;
  static constexpr int STATE_DIM = exp0::STATE_DIM;
  static constexpr int INPUT_DIM = exp0::INPUT_DIM;
  static constexpr Scalar timeStep = 1e-2;
  static constexpr Scalar minRelCost = 1e-3;
  static constexpr Scalar expectedCost = 9.766;
  static constexpr size_t PredictLength = 200;

  using Solver_t = iLQR<Scalar, STATE_DIM, INPUT_DIM, PredictLength>;
  using Problem_t = typename Solver_t::OptimalControlProblem_t;
  using Initializer_t = DefaultInitializer<Scalar, STATE_DIM, INPUT_DIM>;
  using StateVector_t = Vector<Scalar, STATE_DIM>;
  using InputVector_t = Vector<Scalar, INPUT_DIM>;
  using RolloutSettings_t = RolloutSettings<Scalar>;
  using DDPSettings_t = DDPSettings<Scalar>;

  Exp0()
    : problem(exp0::createExp0Problem<Scalar, PredictLength>()),
    initializerPtr(std::make_unique<Initializer_t>())
  {
    setConstantReferenceTrajectory();
  }

  // rollout settings
  RolloutSettings_t rolloutSettings() const {
    RolloutSettings_t rolloutSettings;
    rolloutSettings.timeStep = timeStep;
    return rolloutSettings;
  };

  DDPSettings_t getSettings(SearchStrategyType strategy, bool useFeedbackPolicy = true) const {
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
    testName += "Strategy: " + searchStrategyToString(ddpSettings.strategy_) + " }";
    return testName;
  }

  void performanceIndexTest(const DDPSettings_t& ddpSettings, const PerformanceIndex<Scalar>& performanceIndex) const {
    const auto testName = getTestName(ddpSettings);
    EXPECT_NEAR(performanceIndex.cost, expectedCost, 10.0 * minRelCost) << "MESSAGE: " << testName << ": failed in the total cost test!";
    EXPECT_NEAR(performanceIndex.equalityLagrangian, 0.0, 10.0 * ddpSettings.constraintTolerance_)
      << "MESSAGE: " << testName << ": failed in state-input equality constraint ISE test!";
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
constexpr Exp0::Scalar Exp0::expectedCost;
constexpr size_t Exp0::PredictLength;

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, ddp_feedback_policy) {
  // ddp settings
  auto ddpSettings = getSettings(SearchStrategyType::LINE_SEARCH, true);
  (void)ddpSettings;

  // dynamics and rollout
  exp0::EXP0_Sys1<Scalar> systemDynamics;

  // instantiate
  Solver_t ddp(&systemDynamics, initializerPtr.get());
  exp0::EXP0_Cost<Scalar, static_cast<int>(PredictLength + 1)> cost;
  exp0::EXP0_FinalCost<Scalar, static_cast<int>(PredictLength + 1)> finalCost;
  ddp.optimalControlProblem_.cost.add(cost);
  ddp.optimalControlProblem_.finalCost.add(finalCost);
  ddp.setDesireTrajectory(problem.timeTrajectory, problem.stateTrajectory, problem.inputTrajectory);

  // run ddp
  EXPECT_NO_THROW(ddp.run(startTime, initState));

  // get solution
  const auto& solution = ddp.primalSolution();
  const auto& ctrl = solution.controller_;

  EXPECT_EQ(ctrl.getType(), ControllerType::LINEAR) << "MESSAGE: iLQR solution does not contain a linear feedback policy!";
  EXPECT_DOUBLE_EQ(ctrl.timeStamp_.back(), finalTime) << "MESSAGE: iLQR failed in policy final time of controller!";
  EXPECT_DOUBLE_EQ(solution.timeTrajectory_.back(), finalTime) << "MESSAGE: iLQR failed in policy final time of trajectory!";
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, ddp_feedforward_policy) {
  // ddp settings
  auto ddpSettings = getSettings(ocs2::ddp::Algorithm::SLQ, 2, ocs2::search_strategy::Type::LINE_SEARCH);
  ddpSettings.useFeedbackPolicy_ = false;

  // dynamics and rollout
  ocs2::EXP0_System systemDynamics(referenceManagerPtr);
  ocs2::TimeTriggeredRollout rollout(systemDynamics, rolloutSettings());

  // instantiate
  ocs2::SLQ ddp(ddpSettings, rollout, problem, *initializerPtr);
  ddp.setReferenceManager(referenceManagerPtr);

  // run ddp
  ddp.run(startTime, initState, finalTime);
  // get solution
  const auto solution = ddp.primalSolution(finalTime);
  const auto* ctrlPtr = dynamic_cast<ocs2::FeedforwardController*>(solution.controllerPtr_.get());

  EXPECT_TRUE(ctrlPtr != nullptr) << "MESSAGE: SLQ solution does not contain a feedforward policy!";
  EXPECT_DOUBLE_EQ(ctrlPtr->timeStamp_.back(), finalTime) << "MESSAGE: SLQ failed in policy final time of controller!";
  EXPECT_DOUBLE_EQ(solution.timeTrajectory_.back(), finalTime) << "MESSAGE: SLQ failed in policy final time of trajectory!";
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, ddp_moving_horizon) {
  // ddp settings
  constexpr size_t numThreads = 2;
  const auto ddpSettings = getSettings(ocs2::ddp::Algorithm::SLQ, numThreads, ocs2::search_strategy::Type::LINE_SEARCH);

  // event times
  const ocs2::scalar_array_t eventTimes{ 1.0 };
  const ocs2::size_array_t modeSequence{ 0, 1 };
  referenceManagerPtr = ocs2::getExp0ReferenceManager(eventTimes, modeSequence);

  // dynamics and rollout
  ocs2::EXP0_System systemDynamics(referenceManagerPtr);
  ocs2::TimeTriggeredRollout rollout(systemDynamics, rolloutSettings());

  // instantiate
  ocs2::SLQ ddp(ddpSettings, rollout, problem, *initializerPtr);
  ddp.setReferenceManager(referenceManagerPtr);

  // run SLQ (no active event)
  ocs2::scalar_t startTime = 0.2;
  ocs2::scalar_t finalTime = 0.9;
  EXPECT_NO_THROW(ddp.run(startTime, initState, finalTime));

  // move the time horizon (overlap with one active event)
  startTime = 0.6;
  finalTime = 1.2;
  EXPECT_NO_THROW(ddp.run(startTime, initState, finalTime));

  // move the time horizon (overlap with no active event)
  startTime = 1.1;
  finalTime = 1.5;
  EXPECT_NO_THROW(ddp.run(startTime, initState, finalTime));

  // move the time horizon (no overlap)
  startTime = 1.6;
  finalTime = 2.0;
  EXPECT_NO_THROW(ddp.run(startTime, initState, finalTime));
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_F(Exp0, ddp_hamiltonian) {
  // ddp settings
  constexpr size_t numThreads = 2;
  auto ddpSettings = getSettings(ocs2::ddp::Algorithm::SLQ, numThreads, ocs2::search_strategy::Type::LINE_SEARCH);
  ddpSettings.useFeedbackPolicy_ = true;
  ddpSettings.maxNumIterations_ = 50;
  ddpSettings.minRelCost_ = 1e-9;  // to allow more iterations that the effect of final linesearch is negligible

  // dynamics and rollout
  ocs2::EXP0_System systemDynamics(referenceManagerPtr);
  ocs2::TimeTriggeredRollout rollout(systemDynamics, rolloutSettings());

  // instantiate
  ocs2::SLQ ddp(ddpSettings, rollout, problem, *initializerPtr);
  ddp.setReferenceManager(referenceManagerPtr);

  // run ddp
  ddp.run(startTime, initState, finalTime);
  // get solution
  const auto solution = ddp.primalSolution(finalTime);

  // define precision for tests
  constexpr ocs2::scalar_t precision = 1e-3;

  // get Hamiltonian at current solution
  // expected outcome: true, because the current solution should be optimal
  ocs2::scalar_t time = solution.timeTrajectory_.front();
  ocs2::vector_t state = solution.stateTrajectory_.front();
  ocs2::vector_t input = solution.controllerPtr_->computeInput(time, state);
  auto hamiltonian = ddp.getHamiltonian(time, state, input);
  const ocs2::vector_t dHdu1a = hamiltonian.dfdu;
  EXPECT_TRUE(dHdu1a.isZero(precision)) << "MESSAGE for test 1a: Derivative of Hamiltonian w.r.t. to u is not zero: " << dHdu1a.transpose();

  // evaluate Hamiltonian at different state (but using feedback policy)
  // expected outcome: true, because for a linear system the LQ approximation of H is exact and the linear feedback policy is globally
  // optimal
  ocs2::scalar_t querryTime = solution.timeTrajectory_.front();
  ocs2::vector_t querryState = ocs2::vector_t::Random(solution.stateTrajectory_.front().size());
  ocs2::vector_t querryInput = solution.controllerPtr_->computeInput(querryTime, querryState);
  const ocs2::vector_t dHdu1b = hamiltonian.dfdux * (querryState - state) + hamiltonian.dfduu * (querryInput - input) + hamiltonian.dfdu;
  EXPECT_TRUE(dHdu1b.isZero(precision)) << "MESSAGE for test 1b: Derivative of Hamiltonian w.r.t. to u is not zero: " << dHdu1b.transpose();

  // evaluate Hamiltonian at different input
  // expected outcome: false, because for a linear system the LQ approximation of H is exact and a random input is not optimal
  querryTime = solution.timeTrajectory_.front();
  querryState = solution.stateTrajectory_.front();
  querryInput = ocs2::vector_t::Random(solution.inputTrajectory_.front().size());
  const ocs2::vector_t dHdu1c = hamiltonian.dfdux * (querryState - state) + hamiltonian.dfduu * (querryInput - input) + hamiltonian.dfdu;
  EXPECT_FALSE(dHdu1c.isZero(precision)) << "MESSAGE for test 1c: Derivative of Hamiltonian w.r.t. to u is zero: " << dHdu1c.transpose();

  // get Hamiltonian at different state (but using feedback policy)
  // expected outcome: true, because for a linear system the linear feedback policy is globally optimal
  time = solution.timeTrajectory_.front();
  state = ocs2::vector_t::Random(solution.stateTrajectory_.front().size());
  input = solution.controllerPtr_->computeInput(time, state);
  hamiltonian = ddp.getHamiltonian(time, state, input);
  const ocs2::vector_t dHdu2 = hamiltonian.dfdu;
  EXPECT_TRUE(dHdu2.isZero(precision)) << "MESSAGE for test 2: Derivative of Hamiltonian w.r.t. to u is not zero: " << dHdu2.transpose();

  // get Hamiltonian at different input
  // expected outcome: false, because a random input is not optimal
  time = solution.timeTrajectory_.front();
  state = solution.stateTrajectory_.front();
  input = ocs2::vector_t::Random(solution.inputTrajectory_.front().size());
  hamiltonian = ddp.getHamiltonian(time, state, input);
  const ocs2::vector_t dHdu3 = hamiltonian.dfdu;
  EXPECT_FALSE(dHdu3.isZero(precision)) << "MESSAGE for test 3: Derivative of Hamiltonian w.r.t. to u is zero: " << dHdu3.transpose();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/* Add parameterized test suite */
class Exp0Param : public Exp0, public testing::WithParamInterface<std::tuple<ocs2::search_strategy::Type, size_t>> {
protected:
  ocs2::search_strategy::Type getSearchStrategy() { return std::get<0>(GetParam()); }

  size_t getNumThreads() { return std::get<1>(GetParam()); }
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_P(Exp0Param, SLQ) {
  // ddp settings
  const auto ddpSettings = getSettings(ocs2::ddp::Algorithm::SLQ, getNumThreads(), getSearchStrategy());

  // dynamics and rollout
  ocs2::EXP0_System systemDynamics(referenceManagerPtr);
  ocs2::TimeTriggeredRollout rollout(systemDynamics, rolloutSettings());

  // instantiate
  ocs2::SLQ ddp(ddpSettings, rollout, problem, *initializerPtr);
  ddp.setReferenceManager(referenceManagerPtr);

  if (ddpSettings.displayInfo_ || ddpSettings.displayShortSummary_) {
    std::cerr << "\n" << getTestName(ddpSettings) << "\n";
  }

  // run ddp
  ddp.run(startTime, initState, finalTime);
  // get performance index
  const auto performanceIndex = ddp.getPerformanceIndeces();

  // performanceIndeces test
  performanceIndexTest(ddpSettings, performanceIndex);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
TEST_P(Exp0Param, ILQR) {
  // ddp settings
  const auto ddpSettings = getSettings(ocs2::ddp::Algorithm::ILQR, getNumThreads(), getSearchStrategy());

  // dynamics and rollout
  ocs2::EXP0_System systemDynamics(referenceManagerPtr);
  ocs2::TimeTriggeredRollout rollout(systemDynamics, rolloutSettings());

  // instantiate
  ocs2::ILQR ddp(ddpSettings, rollout, problem, *initializerPtr);
  ddp.setReferenceManager(referenceManagerPtr);

  if (ddpSettings.displayInfo_ || ddpSettings.displayShortSummary_) {
    std::cerr << "\n" << getTestName(ddpSettings) << "\n";
  }

  // run ddp
  ddp.run(startTime, initState, finalTime);
  // get performance index
  const auto performanceIndex = ddp.getPerformanceIndeces();

  // performanceIndeces test
  performanceIndexTest(ddpSettings, performanceIndex);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
INSTANTIATE_TEST_CASE_P(Exp0ParamCase, Exp0Param,
  testing::Combine(testing::ValuesIn({ ocs2::search_strategy::Type::LINE_SEARCH,
                                      ocs2::search_strategy::Type::LEVENBERG_MARQUARDT }),
    testing::ValuesIn({ size_t(1), size_t(3) })), /* num threads */
  [](const testing::TestParamInfo<Exp0Param::ParamType>& info) {
    /* returns test name for gtest summary */
    std::string name;
    name += ocs2::search_strategy::toString(std::get<0>(info.param)) + "__";
    name += std::get<1>(info.param) == 1 ? "SINGLE_THREAD" : "MULTI_THREAD";
    return name;
  });
