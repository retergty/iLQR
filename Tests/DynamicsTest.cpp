/**
 * @file DynamicsTest.cpp
 * @brief 动力学模块测试：LinearSystemDynamics 的 flow map 与线性近似。
 */
#include <gtest/gtest.h>

#include "Controller/LinearController.hpp"
#include "Dynamics/LinearSystemDynamics.hpp"
#include "Matrix/Types.hpp"

// 验证线性系统的流映射返回 Ax + Bu。
TEST(DynamicsTest, LinearSystem_FlowMap) {
  Matrix<double, 2, 2> A = Matrix<double, 2, 2>::Identity();
  Matrix<double, 2, 2> B = Matrix<double, 2, 2>::Identity();
  LinearSystemDynamics<double, 2, 2> dyn(A, B);

  const Vector<double, 2> x{1.0, 2.0};
  const Vector<double, 2> u{0.5, -0.5};

  const Vector<double, 2> xdot = dyn.computeFlowMap(0.0, x, u);
  const Vector<double, 2> expected{1.0 + 0.5, 2.0 - 0.5};
  EXPECT_TRUE(xdot.isApprox(expected, 1e-10));
}

// 验证控制输入为零时流映射退化为 Ax。
TEST(DynamicsTest, LinearSystem_FlowMapZeroInput) {
  Matrix<double, 3, 3> A{{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, -1.0}};
  Matrix<double, 3, 1> B{0.0, 0.0, 1.0};
  LinearSystemDynamics<double, 3, 1> dyn(A, B);

  const Vector<double, 3> x{1.0, 1.0, 1.0};
  const Vector<double, 1> u{0.0};

  const Vector<double, 3> xdot = dyn.computeFlowMap(0.0, x, u);
  EXPECT_DOUBLE_EQ(xdot(0), 1.0);
  EXPECT_DOUBLE_EQ(xdot(1), 2.0);
  EXPECT_DOUBLE_EQ(xdot(2), -1.0);
}

// 验证线性近似返回的 A、B 和名义流值正确。
TEST(DynamicsTest, LinearSystem_LinearApproximation) {
  Matrix<double, 2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
  Matrix<double, 2, 2> B = Matrix<double, 2, 2>::Identity();
  LinearSystemDynamics<double, 2, 2> dyn(A, B);

  const Vector<double, 2> x{1.0, 0.0};
  const Vector<double, 2> u{0.0, 1.0};

  auto approx = dyn.linearApproximation(0.0, x, u);

  EXPECT_TRUE(approx.dfdx.isApprox(A, 1e-10));
  EXPECT_TRUE(approx.dfdu.isApprox(B, 1e-10));
  const Vector<double, 2> expectedF = A * x + B * u;
  EXPECT_TRUE(approx.f.isApprox(expectedF, 1e-10));
}

// 验证受控流映射使用控制器插值得到的输入。
TEST(DynamicsTest, ControlledFlowMap_UsesInterpolatedControllerInput) {
  Matrix<double, 2, 2> A = Matrix<double, 2, 2>::Zero();
  Matrix<double, 2, 2> B = Matrix<double, 2, 2>::Identity();
  LinearSystemDynamics<double, 2, 2> dyn(A, B);

  LinearController<double, 2, 2, 2> controller;
  controller.timeStamp_[0] = 0.0;
  controller.timeStamp_[1] = 1.0;
  controller.biasArray_[0] = {1.0, 2.0};
  controller.biasArray_[1] = {3.0, 4.0};
  controller.gainArray_[0] = Matrix<double, 2, 2>::Zero();
  controller.gainArray_[1] = Matrix<double, 2, 2>::Zero();
  dyn.setController(&controller);

  const Vector<double, 2> x{10.0, -10.0};

  const Vector<double, 2> expected{2.0, 3.0};

  EXPECT_TRUE(dyn.computeFlowMap(0.5, x).isApprox(expected, 1e-10));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
