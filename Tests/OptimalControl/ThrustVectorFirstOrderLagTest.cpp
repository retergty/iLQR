#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "ExampleModels/ThrustVectorFirstOrderLag.hpp"
#include "iLQR/iLQR.hpp"

namespace {
using Scalar = double;
constexpr size_t PredictLength = 15;
constexpr Scalar TimeStep = 0.01;
constexpr Scalar LagAlpha = 0.5;

using Descriptor = iLQRDescriptor<
    Scalar,
    TranscriptionConfig<Dimensions<thrust_vector_first_order_lag::STATE_DIM,
                                   thrust_vector_first_order_lag::INPUT_DIM>,
                        Horizon<PredictLength>, DiscreteDynamics>>;
using Solver = iLQR<Descriptor>;
using ThrustVectorILQR =
    thrust_vector_first_order_lag::ThrustVectorILQR<Scalar, PredictLength>;
using StateVector = typename Solver::StateVector_t;
using InputVector = typename Solver::InputVector_t;

void runOneMpcCycle(ThrustVectorILQR& ilqr,
                    const StateVector& velocityReference, Scalar& currentTime,
                    StateVector& currentState) {
  auto& solver = ilqr.solver();
  auto& problem = ilqr.problem();
  const Vector<Scalar, 3> velocitySetpoint =
      velocityReference.template head<3>();
  const Vector<Scalar, 3> currentVelocity = currentState.template head<3>();
  const Vector<Scalar, 3> commandAccelerationReference =
      currentState.template segment<3>(6);

  ilqr.setDesireTrajectory(currentTime, velocitySetpoint, currentVelocity,
                           commandAccelerationReference);

  ASSERT_NO_THROW(solver.run(currentTime, currentState));
  const InputVector firstInput = solver.primalSolution().inputTrajectory_[0];
  const Vector<Scalar, 3> previousEffectiveAcceleration =
      currentState.template segment<3>(3);
  const Vector<Scalar, 3> previousCommandAcceleration =
      currentState.template segment<3>(6);

  EXPECT_TRUE(firstInput.isAllFinite());
  currentState = problem.dynamicsPtr->computeMap(currentTime, currentState,
                                                 firstInput, TimeStep);
  EXPECT_TRUE(currentState.isAllFinite());

  for (int i = 0; i < 3; ++i) {
    const Scalar commandAcceleration =
        previousCommandAcceleration(i) + firstInput(i);
    const Scalar effectiveAcceleration =
        (Scalar(1) - LagAlpha) * previousEffectiveAcceleration(i) +
        LagAlpha * commandAcceleration;
    EXPECT_NEAR(currentState(i + 3), effectiveAcceleration, Scalar(1e-10));
    EXPECT_NEAR(currentState(i + 6), commandAcceleration, Scalar(1e-10));
  }

  currentTime += TimeStep;
}

Scalar velocityError(const StateVector& state,
                     const StateVector& velocityReference) {
  return (state.template head<3>() - velocityReference.template head<3>())
      .norm();
}

thrust_vector_first_order_lag::ThrustVectorILQRSettings<Scalar>
makeILQRSettings(const size_t maxNumIterations, const Scalar minRelCost) {
  thrust_vector_first_order_lag::ThrustVectorILQRSettings<Scalar> settings;
  settings.ddpSettings.timeStep = TimeStep;
  settings.ddpSettings.maxNumIterations = maxNumIterations;
  settings.ddpSettings.minRelCost = minRelCost;

  settings.ocpSettings.Q.setZero();
  settings.ocpSettings.Q.template topLeftCorner<3, 3>().setIdentity();
  settings.ocpSettings.R =
      Scalar(1e-3) *
      Matrix<Scalar, thrust_vector_first_order_lag::INPUT_DIM,
             thrust_vector_first_order_lag::INPUT_DIM>::Identity();
  settings.ocpSettings.Ra =
      Scalar(1e-4) *
      Matrix<Scalar, thrust_vector_first_order_lag::INPUT_DIM,
             thrust_vector_first_order_lag::INPUT_DIM>::Identity();
  settings.ocpSettings.Qf.setZero();
  settings.ocpSettings.Qf.template topLeftCorner<3, 3>() =
      Scalar(10.0) * Matrix<Scalar, 3, 3>::Identity();
  settings.ocpSettings.weight = Scalar(1);
  settings.ocpSettings.alpha = LagAlpha;
  return settings;
}
}  // namespace

TEST(ThrustVectorFirstOrderLagDynamicsTest, ZeroInputKeepsStatesAtZero) {
  thrust_vector_first_order_lag::ThrustVectorDynamicSystem<Scalar> dynamics(
      LagAlpha);

  StateVector state;
  state.setZero();
  InputVector input;
  input.setZero();

  const StateVector nextState =
      dynamics.computeMap(Scalar(0.0), state, input, TimeStep);

  EXPECT_TRUE(nextState.template head<3>().isZero(Scalar(1e-12)));
  EXPECT_TRUE(nextState.template segment<3>(3).isZero(Scalar(1e-12)));
  EXPECT_TRUE(nextState.template segment<3>(6).isZero(Scalar(1e-12)));
}

TEST(ThrustVectorFirstOrderLagDynamicsTest, InputIsAppliedThroughLag) {
  thrust_vector_first_order_lag::ThrustVectorDynamicSystem<Scalar> dynamics(
      LagAlpha);

  StateVector state;
  state.setZero();
  state(3) = Scalar(0.2);
  state(6) = Scalar(0.6);

  InputVector input;
  input.setZero();
  input(0) = Scalar(0.4);

  const StateVector nextState =
      dynamics.computeMap(Scalar(0.0), state, input, TimeStep);
  const Scalar commandAcceleration = state(6) + input(0);
  const Scalar effectiveAcceleration =
      (Scalar(1) - LagAlpha) * state(3) + LagAlpha * commandAcceleration;

  EXPECT_NEAR(nextState(0), TimeStep * effectiveAcceleration, Scalar(1e-12));
  EXPECT_NEAR(nextState(3), effectiveAcceleration, Scalar(1e-12));
  EXPECT_NEAR(nextState(6), commandAcceleration, Scalar(1e-12));
}

TEST(ThrustVectorFirstOrderLagTrajectoryTest,
     DeploymentInterfaceFiltersVelocityAndWritesCommandReference) {
  auto ilqrSettings = makeILQRSettings(25, Scalar(1e-6));
  ilqrSettings.ocpSettings.referenceTrajectoryAlpha =
      Vector<Scalar, 3>{Scalar(0.5), Scalar(0.25), Scalar(1)};
  auto ilqr =
      thrust_vector_first_order_lag::createThrustVectorFirstOrderLagProblem<
          Scalar, PredictLength>(ilqrSettings);

  const Vector<Scalar, 3> currentVelocity{Scalar(0), Scalar(1), Scalar(-1)};
  const Vector<Scalar, 3> velocitySetpoint{Scalar(2), Scalar(5), Scalar(3)};
  const Vector<Scalar, 3> commandAccelerationReference{
      Scalar(0.1), Scalar(-0.2), Scalar(0.3)};

  ilqr.setDesireTrajectory(Scalar(0), velocitySetpoint, currentVelocity,
                           commandAccelerationReference);

  const auto& targetTrajectory = ilqr.solver().targetTrajectory();
  EXPECT_NEAR(targetTrajectory.stateTrajectory[0](0), Scalar(0.0),
              Scalar(1e-12));
  EXPECT_NEAR(targetTrajectory.stateTrajectory[0](1), Scalar(1.0),
              Scalar(1e-12));
  EXPECT_NEAR(targetTrajectory.stateTrajectory[0](2), Scalar(-1.0),
              Scalar(1e-12));
  EXPECT_NEAR(targetTrajectory.stateTrajectory[1](0), Scalar(1.0),
              Scalar(1e-12));
  EXPECT_NEAR(targetTrajectory.stateTrajectory[1](1), Scalar(2.0),
              Scalar(1e-12));
  EXPECT_NEAR(targetTrajectory.stateTrajectory[1](2), Scalar(3.0),
              Scalar(1e-12));

  for (size_t i = 0; i < PredictLength + 1; ++i) {
    EXPECT_NEAR(targetTrajectory.stateTrajectory[i](6), Scalar(0.1),
                Scalar(1e-6));
    EXPECT_NEAR(targetTrajectory.stateTrajectory[i](7), Scalar(-0.2),
                Scalar(1e-6));
    EXPECT_NEAR(targetTrajectory.stateTrajectory[i](8), Scalar(0.3),
                Scalar(1e-6));
    EXPECT_TRUE(targetTrajectory.inputTrajectory[i].isZero(Scalar(1e-12)));
  }
}

TEST(ThrustVectorFirstOrderLagCostTest,
     TrackCostIncludesCommandAccelerationReference) {
  using TrackCost =
      thrust_vector_first_order_lag::ThrustVectorTrackCost<Scalar, 2>;

  Matrix<Scalar, thrust_vector_first_order_lag::STATE_DIM,
         thrust_vector_first_order_lag::STATE_DIM>
      Q;
  Q.setZero();
  Matrix<Scalar, thrust_vector_first_order_lag::INPUT_DIM,
         thrust_vector_first_order_lag::INPUT_DIM>
      R;
  R.setZero();
  Matrix<Scalar, thrust_vector_first_order_lag::INPUT_DIM,
         thrust_vector_first_order_lag::INPUT_DIM>
      Ra;
  Ra(0, 0) = Scalar(2.0);
  Ra(0, 1) = Scalar(0.1);
  Ra(0, 2) = Scalar(-0.2);
  Ra(1, 0) = Scalar(0.1);
  Ra(1, 1) = Scalar(3.0);
  Ra(1, 2) = Scalar(0.4);
  Ra(2, 0) = Scalar(-0.2);
  Ra(2, 1) = Scalar(0.4);
  Ra(2, 2) = Scalar(4.0);

  TrackCost cost(Q, R, Ra, 0);

  std::array<Scalar, 2> timeTrajectory{Scalar(0), Scalar(1)};
  std::array<StateVector, 2> stateTrajectory;
  std::array<InputVector, 2> inputTrajectory;
  for (size_t i = 0; i < 2; ++i) {
    stateTrajectory[i].setZero();
    inputTrajectory[i].setZero();
  }
  stateTrajectory[0].template segment<3>(6) =
      Vector<Scalar, 3>{Scalar(0.1), Scalar(-0.4), Scalar(0.6)};

  StateVector state;
  state.setZero();
  state.template segment<3>(6) =
      Vector<Scalar, 3>{Scalar(0.5), Scalar(-0.1), Scalar(0.2)};
  InputVector input{Scalar(0.4), Scalar(0.3), Scalar(-0.2)};

  const Vector<Scalar, 3> commandAccelerationDeviation =
      state.template segment<3>(6) + input -
      stateTrajectory[0].template segment<3>(6);
  const Vector<Scalar, 3> weightedCommandAccelerationDeviation =
      Ra * commandAccelerationDeviation;
  const Scalar expectedValue =
      Scalar(0.5) *
      commandAccelerationDeviation.dot(weightedCommandAccelerationDeviation);

  EXPECT_NEAR(cost.getValue(Scalar(0), state, input, timeTrajectory,
                            stateTrajectory, inputTrajectory),
              expectedValue, Scalar(1e-12));

  auto approximation = ScalarFunctionQuadraticApproximation<
      Scalar, thrust_vector_first_order_lag::STATE_DIM,
      thrust_vector_first_order_lag::INPUT_DIM>::Zero();
  cost.addQuadraticApproximation(Scalar(0), state, input, timeTrajectory,
                                 stateTrajectory, inputTrajectory,
                                 approximation);
  EXPECT_NEAR(approximation.f, expectedValue, Scalar(1e-12));
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(approximation.dfdx(i + 6),
                weightedCommandAccelerationDeviation(i), Scalar(1e-12));
    EXPECT_NEAR(approximation.dfdu(i), weightedCommandAccelerationDeviation(i),
                Scalar(1e-12));
    for (int j = 0; j < 3; ++j) {
      EXPECT_NEAR(approximation.dfdxx(i + 6, j + 6), Ra(i, j), Scalar(1e-12));
      EXPECT_NEAR(approximation.dfduu(i, j), Ra(i, j), Scalar(1e-12));
      EXPECT_NEAR(approximation.dfdux(i, j + 6), Ra(i, j), Scalar(1e-12));
    }
  }
}

TEST(ThrustVectorFirstOrderLagCostTest,
     DiagonalTrackCostIncludesCommandAccelerationReference) {
  using TrackCost =
      thrust_vector_first_order_lag::ThrustVectorDiagonalTrackCost<Scalar, 2>;

  Matrix<Scalar, thrust_vector_first_order_lag::STATE_DIM,
         thrust_vector_first_order_lag::STATE_DIM>
      Q;
  Q.setZero();
  Matrix<Scalar, thrust_vector_first_order_lag::INPUT_DIM,
         thrust_vector_first_order_lag::INPUT_DIM>
      R;
  R.setZero();
  Matrix<Scalar, thrust_vector_first_order_lag::INPUT_DIM,
         thrust_vector_first_order_lag::INPUT_DIM>
      Ra;
  Ra.setZero();
  Ra(0, 0) = Scalar(2.0);
  Ra(1, 1) = Scalar(3.0);
  Ra(2, 2) = Scalar(4.0);

  TrackCost cost(Q, R, Ra, 0);

  std::array<Scalar, 2> timeTrajectory{Scalar(0), Scalar(1)};
  std::array<StateVector, 2> stateTrajectory;
  std::array<InputVector, 2> inputTrajectory;
  for (size_t i = 0; i < 2; ++i) {
    stateTrajectory[i].setZero();
    inputTrajectory[i].setZero();
  }
  stateTrajectory[0].template segment<3>(6) =
      Vector<Scalar, 3>{Scalar(0.1), Scalar(-0.4), Scalar(0.6)};

  StateVector state;
  state.setZero();
  state.template segment<3>(6) =
      Vector<Scalar, 3>{Scalar(0.5), Scalar(-0.1), Scalar(0.2)};
  InputVector input{Scalar(0.4), Scalar(0.3), Scalar(-0.2)};

  const Vector<Scalar, 3> commandAccelerationDeviation =
      state.template segment<3>(6) + input -
      stateTrajectory[0].template segment<3>(6);
  Scalar expectedValue = Scalar(0);
  for (int i = 0; i < 3; ++i) {
    expectedValue += Ra(i, i) * commandAccelerationDeviation(i) *
                     commandAccelerationDeviation(i);
  }
  expectedValue *= Scalar(0.5);

  EXPECT_NEAR(cost.getValue(Scalar(0), state, input, timeTrajectory,
                            stateTrajectory, inputTrajectory),
              expectedValue, Scalar(1e-12));

  auto approximation = ScalarFunctionQuadraticApproximation<
      Scalar, thrust_vector_first_order_lag::STATE_DIM,
      thrust_vector_first_order_lag::INPUT_DIM>::Zero();
  cost.addQuadraticApproximation(Scalar(0), state, input, timeTrajectory,
                                 stateTrajectory, inputTrajectory,
                                 approximation);
  EXPECT_NEAR(approximation.f, expectedValue, Scalar(1e-12));
  for (int i = 0; i < 3; ++i) {
    const Scalar weightedCommandAccelerationDeviation =
        Ra(i, i) * commandAccelerationDeviation(i);
    EXPECT_NEAR(approximation.dfdx(i + 6), weightedCommandAccelerationDeviation,
                Scalar(1e-12));
    EXPECT_NEAR(approximation.dfdu(i), weightedCommandAccelerationDeviation,
                Scalar(1e-12));
    EXPECT_NEAR(approximation.dfdxx(i + 6, i + 6), Ra(i, i), Scalar(1e-12));
    EXPECT_NEAR(approximation.dfduu(i, i), Ra(i, i), Scalar(1e-12));
    EXPECT_NEAR(approximation.dfdux(i, i + 6), Ra(i, i), Scalar(1e-12));
  }
}

TEST(ThrustVectorFirstOrderLagMpcTest,
     RecedingHorizonOptimizationReducesVelocityError) {
  const auto ilqrSettings = makeILQRSettings(25, Scalar(1e-6));
  auto ilqr =
      thrust_vector_first_order_lag::createThrustVectorFirstOrderLagProblem<
          Scalar, PredictLength>(ilqrSettings);
  StateVector currentState;
  currentState.setZero();

  StateVector velocityReference;
  velocityReference.setZero();
  velocityReference(0) = Scalar(1.0);

  const Scalar initialVelocityError =
      velocityError(currentState, velocityReference);

  Scalar currentTime = 0.0;
  for (size_t cycle = 0; cycle < 8; ++cycle) {
    runOneMpcCycle(ilqr, velocityReference, currentTime, currentState);
  }

  const Scalar finalVelocityError =
      velocityError(currentState, velocityReference);

  EXPECT_LT(finalVelocityError, initialVelocityError);
  EXPECT_GT(currentState(0), Scalar(0.0));
}

TEST(ThrustVectorFirstOrderLagMpcTest, TracksAndMaintainsVelocityReference) {
  const auto ilqrSettings = makeILQRSettings(30, Scalar(1e-7));
  auto ilqr =
      thrust_vector_first_order_lag::createThrustVectorFirstOrderLagProblem<
          Scalar, PredictLength>(ilqrSettings);
  StateVector currentState;
  currentState.setZero();

  StateVector velocityReference;
  velocityReference.setZero();
  velocityReference(0) = Scalar(1.0);

  Scalar currentTime = 0.0;
  for (size_t cycle = 0; cycle < 35; ++cycle) {
    runOneMpcCycle(ilqr, velocityReference, currentTime, currentState);
  }

  const Scalar trackedError = velocityError(currentState, velocityReference);
  EXPECT_LT(trackedError, Scalar(0.02));

  Scalar maxHoldError = trackedError;
  for (size_t cycle = 0; cycle < 10; ++cycle) {
    runOneMpcCycle(ilqr, velocityReference, currentTime, currentState);
    maxHoldError =
        std::max(maxHoldError, velocityError(currentState, velocityReference));
  }

  EXPECT_LT(maxHoldError, Scalar(0.02));
  EXPECT_NEAR(currentState(0), velocityReference(0), Scalar(0.02));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
