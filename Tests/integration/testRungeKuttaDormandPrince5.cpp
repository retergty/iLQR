/**
 * @file testRungeKuttaDormandPrince5.cpp
 * @brief Dormand-Prince 固定步长积分测试：单步精度、整段积分和反向积分。
 */
#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "RungeKuttaDormandPrince5.hpp"

namespace {
using Scalar = double;
constexpr int XDim = 1;
using StateVector = Vector<Scalar, XDim>;

class ExponentialOde final : public OdeBase<Scalar, XDim> {
 public:
  StateVector computeFlowMap(Scalar t, const StateVector& x) const override {
    (void)t;
    return x;
  }
};
}  // namespace

TEST(RungeKuttaDormandPrince5Test, StepperMatchesExponentialGrowthForOneStep) {
  ExponentialOde system;
  RungeKuttaDormandPrince5Stepper<Scalar, XDim> stepper;

  StateVector x0;
  x0 << 1.0;

  const Scalar t0 = 0.0;
  const Scalar dt = 0.1;
  const StateVector dxdt0 = system.computeFlowMap(t0, x0);

  StateVector x1;
  StateVector dxdt1;
  stepper.doStep(system, x0, dxdt0, t0, dt, x1, dxdt1);

  EXPECT_NEAR(x1(0), std::exp(dt), 1e-9);
  EXPECT_NEAR(dxdt1(0), x1(0), 1e-12);
}

TEST(RungeKuttaDormandPrince5Test, IntegrateConstMatchesAnalyticSolution) {
  ExponentialOde system;
  RungeKuttaDormandPrince5<Scalar, XDim> integrator;

  constexpr size_t NumSamples = 11;
  std::array<Scalar, NumSamples> timeTrajectory{};
  std::array<StateVector, NumSamples> stateTrajectory{};
  Observer<Scalar, XDim> observer(NumSamples, stateTrajectory.data(),
                                  timeTrajectory.data());

  StateVector x0;
  x0 << 1.0;

  integrator.integrateConst(system, observer, x0, 0.0, 1.0, 0.1);

  ASSERT_EQ(observer.getCount(), static_cast<int>(NumSamples));
  for (size_t i = 0; i < NumSamples; ++i) {
    const Scalar t = 0.1 * static_cast<Scalar>(i);
    EXPECT_NEAR(timeTrajectory[i], t, 1e-12);
    EXPECT_NEAR(stateTrajectory[i](0), std::exp(t), 1e-8);
  }
}

TEST(RungeKuttaDormandPrince5Test, IntegrateConstSupportsBackwardTime) {
  ExponentialOde system;
  RungeKuttaDormandPrince5<Scalar, XDim> integrator;

  constexpr size_t NumSamples = 11;
  std::array<Scalar, NumSamples> forwardTimes{};
  std::array<StateVector, NumSamples> forwardStates{};
  Observer<Scalar, XDim> forwardObserver(NumSamples, forwardStates.data(),
                                         forwardTimes.data());

  StateVector x0;
  x0 << 1.0;

  integrator.integrateConst(system, forwardObserver, x0, 0.0, 1.0, 0.1);

  std::array<Scalar, NumSamples> backwardTimes{};
  std::array<StateVector, NumSamples> backwardStates{};
  Observer<Scalar, XDim> backwardObserver(NumSamples, backwardStates.data(),
                                          backwardTimes.data());
  integrator.integrateConst(system, backwardObserver, forwardStates.back(), 1.0,
                            0.0, -0.1);

  ASSERT_EQ(backwardObserver.getCount(), static_cast<int>(NumSamples));
  EXPECT_NEAR(backwardTimes.front(), 1.0, 1e-12);
  EXPECT_NEAR(backwardTimes.back(), 0.0, 1e-12);
  EXPECT_NEAR(backwardStates.front()(0), forwardStates.back()(0), 1e-12);
  EXPECT_NEAR(backwardStates.back()(0), x0(0), 1e-8);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
