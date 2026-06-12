#include <gtest/gtest.h>

#include <cmath>

#include "ExampleModels/ThrustVectorFirstOrderLagConstrain.hpp"

namespace {
using Scalar = double;
using Constraint =
    thrust_vector_first_order_lag::ThrustCommandAccelerationConstraint<Scalar>;
using OCPSettings =
    thrust_vector_first_order_lag::ThrustVectorOCPSettings<Scalar>;
using StateVector = Vector<Scalar, thrust_vector_first_order_lag::STATE_DIM>;
using InputVector = Vector<Scalar, thrust_vector_first_order_lag::INPUT_DIM>;

Constraint makeConstraint() {
  typename Constraint::Config config;
  config.maxTiltAngleRad = Scalar(0.7853981633974483);
  config.axMax = Scalar(10.0);
  config.ayMax = Scalar(10.0);
  config.azMax = Scalar(10.0);
  config.coneSmoothing = Scalar(1e-6);
  return Constraint(config);
}
}  // namespace

TEST(ThrustVectorFirstOrderLagConstrainTest,
     EllipsoidIsCenteredOnTotalCommandAcceleration) {
  const Constraint constraint = makeConstraint();

  StateVector state;
  state.setZero();
  InputVector input;
  input.setZero();

  const auto value = constraint.getValue(Scalar(0), state, input);
  const auto approximation =
      constraint.getQuadraticApproximation(Scalar(0), state, input);

  EXPECT_GT(value(0), Scalar(9.0));
  EXPECT_NEAR(value(1), Scalar(1.0), Scalar(1e-12));

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(approximation.dfdx(1, i + 6), Scalar(0), Scalar(1e-12));
    EXPECT_NEAR(approximation.dfdu(1, i), Scalar(0), Scalar(1e-12));
  }
}

TEST(ThrustVectorFirstOrderLagConstrainTest,
     EllipsoidDerivativesUseTotalCommandAcceleration) {
  const Constraint constraint = makeConstraint();

  StateVector state;
  state.setZero();
  state(6) = Scalar(1.0);
  state(7) = Scalar(2.0);
  state(8) = Scalar(3.0);

  InputVector input;
  input(0) = Scalar(0.5);
  input(1) = Scalar(-1.0);
  input(2) = Scalar(1.0);

  const Scalar inverseMaxSquared = Scalar(1.0) / Scalar(100.0);
  const Vector<Scalar, 3> commandAcceleration{Scalar(1.5), Scalar(1.0),
                                              Scalar(4.0)};
  const Scalar expectedEllipsoid =
      Scalar(1.0) - inverseMaxSquared * commandAcceleration.squaredNorm();

  const auto value = constraint.getValue(Scalar(0), state, input);
  const auto approximation =
      constraint.getQuadraticApproximation(Scalar(0), state, input);

  EXPECT_NEAR(value(1), expectedEllipsoid, Scalar(1e-12));

  for (int i = 0; i < 3; ++i) {
    const Scalar expectedDerivative =
        -Scalar(2.0) * inverseMaxSquared * commandAcceleration(i);
    EXPECT_NEAR(approximation.dfdx(1, i + 6), expectedDerivative,
                Scalar(1e-12));
    EXPECT_NEAR(approximation.dfdu(1, i), expectedDerivative, Scalar(1e-12));
    EXPECT_NEAR(approximation.dfdxx[1](i + 6, i + 6),
                -Scalar(2.0) * inverseMaxSquared, Scalar(1e-12));
    EXPECT_NEAR(approximation.dfduu[1](i, i), -Scalar(2.0) * inverseMaxSquared,
                Scalar(1e-12));
    EXPECT_NEAR(approximation.dfdux[1](i, i + 6),
                -Scalar(2.0) * inverseMaxSquared, Scalar(1e-12));
  }
}

TEST(ThrustVectorFirstOrderLagConstrainTest,
     OcpSettingsExposeConstraintConfiguration) {
  OCPSettings settings;
  settings.constraintSettings.constraint.axMax = Scalar(5.0);
  settings.constraintSettings.constraint.ayMax = Scalar(6.0);
  settings.constraintSettings.constraint.azMax = Scalar(7.0);
  settings.constraintSettings.zMinConstraint.zMin = Scalar(-2.0);
  settings.constraintSettings.constraintPenalty.scale = Scalar(20.0);
  settings.constraintSettings.zMinPenalty.stepSize = Scalar(0.5);

  const Constraint constraint(settings.constraintSettings.constraint);

  StateVector state;
  state.setZero();
  InputVector input;
  input.setZero();
  input(0) = Scalar(5.0);

  const auto value = constraint.getValue(Scalar(0), state, input);

  EXPECT_NEAR(value(1), Scalar(0.0), Scalar(1e-12));
  EXPECT_NEAR(settings.constraintSettings.zMinConstraint.zMin, Scalar(-2.0),
              Scalar(1e-12));
  EXPECT_NEAR(settings.constraintSettings.constraintPenalty.scale, Scalar(20.0),
              Scalar(1e-12));
  EXPECT_NEAR(settings.constraintSettings.zMinPenalty.stepSize, Scalar(0.5),
              Scalar(1e-12));
}
