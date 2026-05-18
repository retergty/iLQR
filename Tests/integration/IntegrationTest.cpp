/**
 * @file IntegrationTest.cpp
 * @brief Integration 抽象接口测试：通过 IntegratorBase 指针驱动固定步长积分。
 */
#include <array>
#include <memory>

#include <gtest/gtest.h>

#include "Integration.hpp"
#include "LinearController.hpp"
#include "LinearSystemDynamics.hpp"
#include "RungeKuttaDormandPrince5.hpp"

namespace {
using Scalar = double;
constexpr int XDim = 2;
constexpr int UDim = 1;
constexpr size_t ControllerLen = 2;
using StateVector = Vector<Scalar, XDim>;
using InputVector = Vector<Scalar, UDim>;
using StateMatrix = Matrix<Scalar, XDim, XDim>;
using InputMatrix = Matrix<Scalar, XDim, UDim>;
using Controller = LinearController<Scalar, XDim, UDim, ControllerLen>;
using System = LinearSystemDynamics<Scalar, XDim, UDim>;

System makeSecondOrderSystem() {
  StateMatrix A;
  A << -2.0, -1.0, 1.0, 0.0;
  InputMatrix B;
  B << 1.0, 0.0;
  return System(A, B);
}

Controller makeConstantController(Scalar t0, Scalar t1,
                                  const InputVector &input) {
  std::array<Scalar, ControllerLen> controllerTime = {t0, t1};
  std::array<InputVector, ControllerLen> controllerBias = {input, input};
  std::array<Matrix<Scalar, UDim, XDim>, ControllerLen> controllerGain = {
      Matrix<Scalar, UDim, XDim>::Zero(), Matrix<Scalar, UDim, XDim>::Zero()};
  return Controller(controllerTime, controllerBias, controllerGain);
}
} // namespace

TEST(IntegrationTest, IntegratorBaseIntegratesControlledSecondOrderSystem) {
  constexpr Scalar t0 = 0.0;
  constexpr Scalar t1 = 10.0;
  constexpr Scalar dt = 0.05;
  constexpr size_t NumSamples = 201;

  System system = makeSecondOrderSystem();
  InputVector constantInput = InputVector::Ones();
  Controller controller = makeConstantController(t0, t1, constantInput);
  system.setController(&controller);

  std::unique_ptr<IntegratorBase<Scalar, XDim>> integrator =
      std::make_unique<RungeKuttaDormandPrince5<Scalar, XDim>>();

  std::array<Scalar, NumSamples> timeTrajectory{};
  std::array<StateVector, NumSamples> stateTrajectory{};
  Observer<Scalar, XDim> observer(NumSamples, stateTrajectory.data(),
                                  timeTrajectory.data());

  const StateVector x0 = StateVector::Zero();
  integrator->integrateConst(system, observer, x0, t0, t1, dt);

  ASSERT_EQ(observer.getCount(), static_cast<int>(NumSamples));
  EXPECT_NEAR(timeTrajectory.front(), t0, 1e-12);
  EXPECT_NEAR(timeTrajectory.back(), t1, 1e-12);
  EXPECT_TRUE(stateTrajectory.front().isApprox(x0, 1e-12));
  EXPECT_NEAR(stateTrajectory.back()(0), 0.0, 1e-3);
  EXPECT_NEAR(stateTrajectory.back()(1), 1.0, 1e-3);
}

TEST(IntegrationTest, IntegratorBaseSupportsBackwardIntegration) {
  constexpr Scalar t0 = 0.0;
  constexpr Scalar t1 = 1.0;
  constexpr Scalar dt = 0.05;
  constexpr size_t NumSamples = 21;

  System system = makeSecondOrderSystem();
  InputVector constantInput = InputVector::Ones();
  Controller controller = makeConstantController(t0, t1, constantInput);
  system.setController(&controller);

  RungeKuttaDormandPrince5<Scalar, XDim> integrator;
  std::array<Scalar, NumSamples> forwardTimes{};
  std::array<StateVector, NumSamples> forwardStates{};
  Observer<Scalar, XDim> forwardObserver(NumSamples, forwardStates.data(),
                                         forwardTimes.data());

  const StateVector x0 = StateVector::Zero();
  integrator.integrateConst(system, forwardObserver, x0, t0, t1, dt);

  std::array<Scalar, NumSamples> backwardTimes{};
  std::array<StateVector, NumSamples> backwardStates{};
  Observer<Scalar, XDim> backwardObserver(NumSamples, backwardStates.data(),
                                          backwardTimes.data());
  integrator.integrateConst(system, backwardObserver, forwardStates.back(), t1,
                            t0, -dt);

  ASSERT_EQ(backwardObserver.getCount(), static_cast<int>(NumSamples));
  EXPECT_NEAR(backwardTimes.front(), t1, 1e-12);
  EXPECT_NEAR(backwardTimes.back(), t0, 1e-12);
  EXPECT_TRUE(backwardStates.front().isApprox(forwardStates.back(), 1e-12));
  EXPECT_LT((backwardStates.back() - x0).norm(), 2e-3);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
