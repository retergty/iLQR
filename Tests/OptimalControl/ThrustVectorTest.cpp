#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "Initialization/Initializer.hpp"
#include "Models/ThrustVector.hpp"
#include "iLQR/iLQR.hpp"

namespace {
using Scalar = double;
constexpr size_t PredictLength = 15;
constexpr Scalar TimeStep = 0.01;
constexpr Scalar Mass = 1.0;
constexpr Scalar Gravity = 9.8;

using Descriptor = iLQRDescriptor<
    Scalar, TranscriptionConfig<
                Dimensions<thrust_vector::STATE_DIM, thrust_vector::INPUT_DIM>,
                Horizon<PredictLength>, DiscreteDynamics>>;
using Solver = iLQR<Descriptor>;
using StateVector = typename Solver::StateVector_t;
using InputVector = typename Solver::InputVector_t;
using TimeTrajectory = typename Solver::TimeTrajectory_t;
using StateTrajectory = typename Solver::StateTrajectory_t;
using InputTrajectory = typename Solver::InputTrajectory_t;

InputVector hoverInput() {
  return InputVector{Scalar(0.0), Scalar(0.0), Mass * Gravity};
}

class HoverInitializer final
    : public Initializer<Scalar, thrust_vector::STATE_DIM,
                         thrust_vector::INPUT_DIM> {
 public:
  void compute(const Scalar time, const StateVector& state,
               const Scalar nextTime, InputVector& input,
               StateVector& nextState) override {
    (void)time;
    const Scalar dt = nextTime - time;
    input = hoverInput();

    nextState.template head<3>() = state.template head<3>();
    nextState(2) += dt * (input(2) / Mass - Gravity);
    nextState.template tail<3>() = input;
  }
};

void configureVelocityReference(const Scalar initTime,
                                const StateVector& velocityReference,
                                TimeTrajectory& timeTrajectory,
                                StateTrajectory& stateTrajectory,
                                InputTrajectory& inputTrajectory) {
  for (size_t i = 0; i < PredictLength + 1; ++i) {
    timeTrajectory[i] = initTime + static_cast<Scalar>(i) * TimeStep;
    stateTrajectory[i].setZero();
    stateTrajectory[i].template head<3>() =
        velocityReference.template head<3>();
    stateTrajectory[i].template tail<3>() = hoverInput();
    inputTrajectory[i] = hoverInput();
  }
}

void runOneMpcCycle(Solver& solver, Solver::OptimalControlProblem_t& problem,
                    const StateVector& velocityReference, Scalar& currentTime,
                    StateVector& currentState) {
  TimeTrajectory timeTrajectory;
  StateTrajectory stateTrajectory;
  InputTrajectory inputTrajectory;
  configureVelocityReference(currentTime, velocityReference, timeTrajectory,
                             stateTrajectory, inputTrajectory);
  solver.setDesireTrajectory(timeTrajectory, stateTrajectory, inputTrajectory);

  ASSERT_NO_THROW(solver.run(currentTime, currentState));
  const InputVector firstInput = solver.primalSolution().inputTrajectory_[0];

  EXPECT_TRUE(firstInput.isAllFinite());
  currentState = problem.dynamicsPtr->computeMap(currentTime, currentState,
                                                 firstInput, TimeStep);
  EXPECT_TRUE(currentState.isAllFinite());
  const InputVector currentThrust = currentState.template tail<3>();
  EXPECT_TRUE(currentThrust.isApprox(firstInput, 1e-10));

  currentTime += TimeStep;
}

Scalar velocityError(const StateVector& state,
                     const StateVector& velocityReference) {
  return (state.template head<3>() - velocityReference.template head<3>())
      .norm();
}
}  // namespace

TEST(ThrustVectorMpcTest, RecedingHorizonOptimizationReducesVelocityError) {
  auto& problem =
      thrust_vector::createThrustVectorProblem<Scalar, PredictLength>(Mass);
  HoverInitializer initializer;

  DDPSettings<Scalar> settings;
  settings.timeStep_ = TimeStep;
  settings.maxNumIterations_ = 20;
  settings.minRelCost_ = 1e-6;

  Solver solver(settings, problem, &initializer);

  StateVector currentState;
  currentState.setZero();
  currentState.template tail<3>() = hoverInput();

  StateVector velocityReference;
  velocityReference.setZero();
  velocityReference(0) = Scalar(1.0);

  const Scalar initialVelocityError =
      velocityError(currentState, velocityReference);

  Scalar currentTime = 0.0;
  for (size_t cycle = 0; cycle < 5; ++cycle) {
    runOneMpcCycle(solver, problem, velocityReference, currentTime,
                   currentState);
  }

  const Scalar finalVelocityError =
      velocityError(currentState, velocityReference);

  EXPECT_LT(finalVelocityError, initialVelocityError);
  EXPECT_GT(currentState(0), Scalar(0.0));
}

TEST(ThrustVectorMpcTest, TracksAndMaintainsVelocityReference) {
  auto& problem =
      thrust_vector::createThrustVectorProblem<Scalar, PredictLength>(Mass);
  HoverInitializer initializer;

  DDPSettings<Scalar> settings;
  settings.timeStep_ = TimeStep;
  settings.maxNumIterations_ = 25;
  settings.minRelCost_ = 1e-7;

  Solver solver(settings, problem, &initializer);

  StateVector currentState;
  currentState.setZero();
  currentState.template tail<3>() = hoverInput();

  StateVector velocityReference;
  velocityReference.setZero();
  velocityReference(0) = Scalar(1.0);

  Scalar currentTime = 0.0;
  for (size_t cycle = 0; cycle < 25; ++cycle) {
    runOneMpcCycle(solver, problem, velocityReference, currentTime,
                   currentState);
  }

  const Scalar trackedError = velocityError(currentState, velocityReference);
  EXPECT_LT(trackedError, Scalar(0.01));

  Scalar maxHoldError = trackedError;
  for (size_t cycle = 0; cycle < 10; ++cycle) {
    runOneMpcCycle(solver, problem, velocityReference, currentTime,
                   currentState);
    maxHoldError =
        std::max(maxHoldError, velocityError(currentState, velocityReference));
  }

  EXPECT_LT(maxHoldError, Scalar(0.01));
  EXPECT_NEAR(currentState(0), velocityReference(0), Scalar(0.01));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
