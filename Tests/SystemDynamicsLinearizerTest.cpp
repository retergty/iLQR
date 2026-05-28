/**
 * @file SystemDynamicsLinearizerTest.cpp
 * @brief SystemDynamicsLinearizer 测试：验证数值线性化结果与解析 Jacobian
 * 一致。
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>

#include "LinearSystemDynamics.hpp"
#include "SystemDynamicsBase.hpp"
#include "SystemDynamicsLinearizer.hpp"
#include "Types.hpp"

namespace {

// 数值线性化与解析 Jacobian 对比时允许的误差。
constexpr double kTolerance = 1e-5;
// 有限差分基础步长。
constexpr double kEpsilon = 1e-10;
constexpr int kStateDim = 2;
constexpr int kInputDim = 1;
constexpr double kPi = 3.14159265358979323846;

using Scalar = double;
using StateVector = Eigen::Vector<Scalar, kStateDim>;
using InputVector = Eigen::Vector<Scalar, kInputDim>;
using StateMatrix = Eigen::Matrix<Scalar, kStateDim, kStateDim>;
using InputMatrix = Eigen::Matrix<Scalar, kStateDim, kInputDim>;
using LinearApproximation =
    VectorFunctionLinearApproximation<Scalar, kStateDim, kStateDim, kInputDim>;
using SystemDynamics = SystemDynamicsBase<Scalar, kStateDim, kInputDim>;
using LinearSystem = LinearSystemDynamics<Scalar, kStateDim, kInputDim>;

// 比较两个系统在线性化点处的 A、B 是否足够接近。
bool derivativeChecker(SystemDynamics& sys1, SystemDynamics& sys2,
                       Scalar tolerance, Scalar t, const StateVector& x,
                       const InputVector& u);

// 简单摆系统，theta = 0 对应竖直向上平衡点。
class PendulumSystem final : public SystemDynamics {
 public:
  PendulumSystem() = default;
  ~PendulumSystem() override = default;

  StateVector computeFlowMap(Scalar t, const StateVector& x,
                             const InputVector& u) override {
    (void)t;
    StateVector dfdt;
    dfdt << x(1), std::sin(x(0)) + 0.1 * u(0);
    return dfdt;
  }

  LinearApproximation linearApproximation(Scalar t, const StateVector& x,
                                          const InputVector& u) override {
    (void)t;
    LinearApproximation linearDynamics;
    linearDynamics.f = computeFlowMap(t, x, u);
    linearDynamics.dfdx << 0.0, 1.0, std::cos(x(0)), 0.0;
    linearDynamics.dfdu << 0.0, 0.1;
    return linearDynamics;
  }
};

// 验证辅助比较函数能识别 Jacobian 不一致的系统。
TEST(testSystemDynamicsLinearizer, testDerivativeChecker) {
  const Scalar time = 0.0;
  const StateVector state = StateVector::Zero();
  const InputVector input = InputVector::Zero();

  StateMatrix A;
  InputMatrix B;
  A << 0.6, 1.2, -0.8, 3.4;
  B << 1.0, 1.0;

  LinearSystem linSys(A, B);

  A(0, 0) = 0.0;
  B(0, 0) = 0.0;
  LinearSystem alteredSys(A, B);

  ASSERT_FALSE(
      derivativeChecker(linSys, alteredSys, kTolerance, time, state, input));
}

// 验证线性系统经有限差分线性化后与其解析 Jacobian 一致。
TEST(testSystemDynamicsLinearizer, testLinearSystem) {
  const Scalar time = 0.0;
  const StateVector state = StateVector::Zero();
  const InputVector input = InputVector::Zero();

  StateMatrix A;
  InputMatrix B;
  A << 0.6, 1.2, -0.8, 3.4;
  B << 1.0, 1.0;

  LinearSystem linSys(A, B);

  SystemDynamicsLinearizer<Scalar, kStateDim, kInputDim> linearizedSys(
      &linSys,
      /*doubleSidedDerivative=*/true,
      /*isSecondOrderSystem=*/false, kEpsilon);
  ASSERT_TRUE(
      derivativeChecker(linSys, linearizedSys, kTolerance, time, state, input));
}

// 在一组摆角采样点上验证数值线性化与解析摆动力学 Jacobian 一致。
TEST(testSystemDynamicsLinearizer, testPendulum) {
  std::srand(0);

  constexpr std::size_t kDivisions = 1000;
  const Scalar maxDeg = 180.0;
  const Scalar toRads = kPi / 180.0;
  const Scalar t = 0.0;

  std::array<StateVector, kDivisions> testStates{};
  for (std::size_t i = 0; i < kDivisions; ++i) {
    const Scalar alpha =
        static_cast<Scalar>(i) / static_cast<Scalar>(kDivisions - 1);
    testStates[i] << alpha * toRads * maxDeg, 0.0;
  }

  const InputVector input = InputVector::Random();

  PendulumSystem nonLinSys;

  SystemDynamicsLinearizer<Scalar, kStateDim, kInputDim> linearizedSys(
      &nonLinSys,
      /*doubleSidedDerivative=*/true,
      /*isSecondOrderSystem=*/false, kEpsilon);
  for (const auto& state : testStates) {
    ASSERT_TRUE(derivativeChecker(nonLinSys, linearizedSys, kTolerance, t,
                                  state, input));
  }
}

// 仅比较对状态和输入的一阶 Jacobian；名义流值在这里不参与判定。
bool derivativeChecker(SystemDynamics& sys1, SystemDynamics& sys2,
                       Scalar tolerance, Scalar t, const StateVector& x,
                       const InputVector& u) {
  const auto derivatives1 = sys1.linearApproximation(t, x, u);
  const auto derivatives2 = sys2.linearApproximation(t, x, u);
  const Scalar AError =
      (derivatives1.dfdx - derivatives2.dfdx).lpNorm<Eigen::Infinity>();
  const Scalar BError =
      (derivatives1.dfdu - derivatives2.dfdu).lpNorm<Eigen::Infinity>();
  return tolerance > std::max(AError, BError);
}

}  // namespace