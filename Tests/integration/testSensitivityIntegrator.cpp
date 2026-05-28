/**
 * @file testSensitivityIntegrator.cpp
 * @brief 动力学离散器测试：验证 Euler、RK2、RK4 的离散化与一阶敏感度实现。
 */
#include <gtest/gtest.h>

#include <array>

#include "LinearController.hpp"
#include "LinearSystemDynamics.hpp"
#include "RungeKuttaDormandPrince5.hpp"
#include "SensitivityIntegrator.hpp"

namespace {
using Scalar = double;
constexpr int XDim = 2;
constexpr int UDim = 1;
constexpr size_t ControllerLen = 2;

using StateVector = Vector<Scalar, XDim>;
using InputVector = Vector<Scalar, UDim>;
using StateMatrix = Matrix<Scalar, XDim, XDim>;
using InputMatrix = Matrix<Scalar, XDim, UDim>;
using LinearApproximation =
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>;
using Controller = LinearController<Scalar, XDim, UDim, ControllerLen>;
using System = LinearSystemDynamics<Scalar, XDim, UDim>;

System makeSystem() {
  StateMatrix A{{-2.0, -1.0}, {1.0, 0.0}};
  InputMatrix B{1.0, 0.0};
  return System(A, B);
}

Controller makeConstantController(Scalar t0, Scalar t1,
                                  const InputVector& input) {
  std::array<Scalar, ControllerLen> controllerTime = {t0, t1};
  std::array<InputVector, ControllerLen> controllerBias = {input, input};
  std::array<Matrix<Scalar, UDim, XDim>, ControllerLen> controllerGain = {
      Matrix<Scalar, UDim, XDim>::Zero(), Matrix<Scalar, UDim, XDim>::Zero()};
  return Controller(controllerTime, controllerBias, controllerGain);
}
}  // namespace

TEST(SensitivityIntegratorTest, EulerSensitivityMatchesManualConstruction) {
  System system = makeSystem();
  EulerDynamicsDiscretizer<Scalar, XDim, UDim> discretizer;

  const Scalar t = 0.5;
  const Scalar dt = 0.1;
  const StateVector x{0.2, -0.4};
  const InputVector u{0.7};

  const LinearApproximation k1 = system.linearApproximation(t, x, u);

  LinearApproximation expected;
  expected.dfdx = StateMatrix::Identity() + dt * k1.dfdx;
  expected.dfdu = dt * k1.dfdu;
  expected.f = x + dt * k1.f;

  const auto actualDynamics = discretizer.discretize(system, t, x, u, dt);
  const auto actualApproximation =
      discretizer.sensitivityDiscretize(system, t, x, u, dt);

  EXPECT_TRUE(actualDynamics.isApprox(expected.f, 1e-12));
  EXPECT_TRUE(actualApproximation.f.isApprox(expected.f, 1e-12));
  EXPECT_TRUE(actualApproximation.dfdx.isApprox(expected.dfdx, 1e-12));
  EXPECT_TRUE(actualApproximation.dfdu.isApprox(expected.dfdu, 1e-12));
}

TEST(SensitivityIntegratorTest, RK2SensitivityMatchesManualConstruction) {
  System system = makeSystem();
  EK2DynamicsDiscretizer<Scalar, XDim, UDim> discretizer;

  const Scalar t = 0.5;
  const Scalar dt = 0.1;
  const Scalar dtHalve = dt / 2.0;
  const StateVector x{0.2, -0.4};
  const InputVector u{0.7};

  const LinearApproximation k1 = system.linearApproximation(t, x, u);
  const LinearApproximation k2 =
      system.linearApproximation(t + dt, x + dt * k1.f, u);

  LinearApproximation expected;
  expected.dfdx = StateMatrix::Identity() + dtHalve * k1.dfdx +
                  dtHalve * (k2.dfdx + dt * k2.dfdx * k1.dfdx);
  expected.dfdu =
      dtHalve * k1.dfdu + dtHalve * (k2.dfdu + dt * k2.dfdx * k1.dfdu);
  expected.f = x + dtHalve * k1.f + dtHalve * k2.f;

  const auto actualDynamics = discretizer.discretize(system, t, x, u, dt);
  const auto actualApproximation =
      discretizer.sensitivityDiscretize(system, t, x, u, dt);

  EXPECT_TRUE(actualDynamics.isApprox(expected.f, 1e-12));
  EXPECT_TRUE(actualApproximation.f.isApprox(expected.f, 1e-12));
  EXPECT_TRUE(actualApproximation.dfdx.isApprox(expected.dfdx, 1e-12));
  EXPECT_TRUE(actualApproximation.dfdu.isApprox(expected.dfdu, 1e-12));
}

TEST(SensitivityIntegratorTest, RK4SensitivityMatchesManualConstruction) {
  System system = makeSystem();
  EK4DynamicsDiscretizer<Scalar, XDim, UDim> discretizer;

  const Scalar t = 0.5;
  const Scalar dt = 0.1;
  const Scalar dtHalve = dt / 2.0;
  const Scalar dtThird = dt / 3.0;
  const Scalar dtSixth = dt / 6.0;
  const StateVector x{0.2, -0.4};
  const InputVector u{0.7};

  const LinearApproximation k1 = system.linearApproximation(t, x, u);
  const LinearApproximation k2 =
      system.linearApproximation(t + dtHalve, x + dtHalve * k1.f, u);
  const LinearApproximation k3 =
      system.linearApproximation(t + dtHalve, x + dtHalve * k2.f, u);
  const LinearApproximation k4 =
      system.linearApproximation(t + dt, x + dt * k3.f, u);

  const StateMatrix dk1dxk = k1.dfdx;
  const StateMatrix dk2dxk = k2.dfdx + dtHalve * k2.dfdx * dk1dxk;
  const StateMatrix dk3dxk = k3.dfdx + dtHalve * k3.dfdx * dk2dxk;
  const StateMatrix dk4dxk = k4.dfdx + dt * k4.dfdx * dk3dxk;

  const InputMatrix dk1duk = k1.dfdu;
  const InputMatrix dk2duk = k2.dfdu + dtHalve * k2.dfdx * dk1duk;
  const InputMatrix dk3duk = k3.dfdu + dtHalve * k3.dfdx * dk2duk;
  const InputMatrix dk4duk = k4.dfdu + dt * k4.dfdx * dk3duk;

  LinearApproximation expected;
  expected.dfdx = StateMatrix::Identity() + dtSixth * dk1dxk +
                  dtThird * dk2dxk + dtThird * dk3dxk + dtSixth * dk4dxk;
  expected.dfdu =
      dtSixth * dk1duk + dtThird * dk2duk + dtThird * dk3duk + dtSixth * dk4duk;
  expected.f =
      x + dtSixth * k1.f + dtThird * k2.f + dtThird * k3.f + dtSixth * k4.f;

  const auto actualDynamics = discretizer.discretize(system, t, x, u, dt);
  const auto actualApproximation =
      discretizer.sensitivityDiscretize(system, t, x, u, dt);

  EXPECT_TRUE(actualDynamics.isApprox(expected.f, 1e-12));
  EXPECT_TRUE(actualApproximation.f.isApprox(expected.f, 1e-12));
  EXPECT_TRUE(actualApproximation.dfdx.isApprox(expected.dfdx, 1e-12));
  EXPECT_TRUE(actualApproximation.dfdu.isApprox(expected.dfdu, 1e-12));
}

TEST(SensitivityIntegratorTest,
     RK4DiscretizationMatchesDormandPrinceForConstantInput) {
  System system = makeSystem();
  EK4DynamicsDiscretizer<Scalar, XDim, UDim> discretizer;
  RungeKuttaDormandPrince5<Scalar, XDim> integrator;

  const Scalar t = 0.5;
  const Scalar dt = 0.01;
  const StateVector x{0.2, -0.4};
  const InputVector u{0.7};

  Controller controller = makeConstantController(t, t + dt, u);
  system.setController(&controller);

  std::array<Scalar, 2> timeTrajectory{};
  std::array<StateVector, 2> stateTrajectory{};
  Observer<Scalar, XDim> observer(2, stateTrajectory.data(),
                                  timeTrajectory.data());
  integrator.integrateConst(system, observer, x, t, t + dt, dt);

  const StateVector rk4ForwardDynamics =
      discretizer.discretize(system, t, x, u, dt);

  ASSERT_EQ(observer.getCount(), 2);
  EXPECT_NEAR(timeTrajectory.front(), t, 1e-12);
  EXPECT_NEAR(timeTrajectory.back(), t + dt, 1e-12);
  EXPECT_TRUE(stateTrajectory.back().isApprox(rk4ForwardDynamics, 1e-9));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
