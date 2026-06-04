#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "ExampleModels/ThrustVector.hpp"
#include "Initialization/Initializer.hpp"
#include "iLQR/iLQR.hpp"

namespace {
using Scalar = double;
constexpr size_t PredictLength = 15;
constexpr Scalar TimeStep = 0.01;

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
  return InputVector{Scalar(0.0), Scalar(0.0), Scalar(0.0)};
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

    nextState.template head<3>() = state.template head<3>() + dt * input;
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
  const InputVector currentAcceleration = currentState.template tail<3>();
  EXPECT_TRUE(currentAcceleration.isApprox(firstInput, 1e-10));

  currentTime += TimeStep;
}

Scalar velocityError(const StateVector& state,
                     const StateVector& velocityReference) {
  return (state.template head<3>() - velocityReference.template head<3>())
      .norm();
}
}  // namespace

TEST(ThrustVectorDynamicsTest, ZeroInputMaintainsVelocity) {
  thrust_vector::ThrustVectorDynamicSystem<Scalar> dynamics;

  StateVector state;
  state.setZero();
  InputVector input;
  input.setZero();

  const StateVector nextState =
      dynamics.computeMap(Scalar(0.0), state, input, TimeStep);

  EXPECT_TRUE(nextState.template head<3>().isZero(Scalar(1e-12)));
  EXPECT_NEAR(nextState(3), input(0), Scalar(1e-12));
  EXPECT_NEAR(nextState(4), input(1), Scalar(1e-12));
  EXPECT_NEAR(nextState(5), input(2), Scalar(1e-12));
}

TEST(ThrustVectorMpcTest, RecedingHorizonOptimizationReducesVelocityError) {
  auto& problem =
      thrust_vector::createThrustVectorProblem<Scalar, PredictLength>();
  HoverInitializer initializer;

  DDPSettings<Scalar> settings;
  settings.timeStep = TimeStep;
  settings.maxNumIterations = 20;
  settings.minRelCost = 1e-6;

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
      thrust_vector::createThrustVectorProblem<Scalar, PredictLength>();
  HoverInitializer initializer;

  DDPSettings<Scalar> settings;
  settings.timeStep = TimeStep;
  settings.maxNumIterations = 25;
  settings.minRelCost = 1e-7;

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
