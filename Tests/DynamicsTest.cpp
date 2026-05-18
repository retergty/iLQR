/**
 * @file DynamicsTest.cpp
 * @brief 动力学模块测试：LinearSystemDynamics 的 flow map 与线性近似。
 */
#include <gtest/gtest.h>

#include "LinearController.hpp"
#include "LinearSystemDynamics.hpp"
#include "Types.hpp"

// 验证线性系统的流映射返回 Ax + Bu。
TEST(DynamicsTest, LinearSystem_FlowMap) {
  Eigen::Matrix2d A = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d B = Eigen::Matrix2d::Identity();
  LinearSystemDynamics<double, 2, 2> dyn(A, B);

  Eigen::Vector2d x, u;
  x << 1.0, 2.0;
  u << 0.5, -0.5;

  Eigen::Vector2d xdot = dyn.computeFlowMap(0.0, x, u);
  Eigen::Vector2d expected;
  expected << 1.0 + 0.5, 2.0 - 0.5;
  EXPECT_TRUE(xdot.isApprox(expected, 1e-10));
}

// 验证控制输入为零时流映射退化为 Ax。
TEST(DynamicsTest, LinearSystem_FlowMapZeroInput) {
  Eigen::Matrix3d A;
  A << 1, 0, 0, 0, 2, 0, 0, 0, -1;
  Eigen::Matrix<double, 3, 1> B;
  B << 0, 0, 1;
  LinearSystemDynamics<double, 3, 1> dyn(A, B);

  Eigen::Vector3d x;
  x << 1.0, 1.0, 1.0;
  Eigen::Vector<double, 1> u;
  u << 0.0;

  Eigen::Vector3d xdot = dyn.computeFlowMap(0.0, x, u);
  EXPECT_DOUBLE_EQ(xdot(0), 1.0);
  EXPECT_DOUBLE_EQ(xdot(1), 2.0);
  EXPECT_DOUBLE_EQ(xdot(2), -1.0);
}

// 验证线性近似返回的 A、B 和名义流值正确。
TEST(DynamicsTest, LinearSystem_LinearApproximation) {
  Eigen::Matrix2d A, B;
  A << 1, 2, 3, 4;
  B << 1, 0, 0, 1;
  LinearSystemDynamics<double, 2, 2> dyn(A, B);

  Eigen::Vector2d x, u;
  x << 1.0, 0.0;
  u << 0.0, 1.0;

  auto approx = dyn.linearApproximation(0.0, x, u);

  EXPECT_TRUE(approx.dfdx.isApprox(A, 1e-10));
  EXPECT_TRUE(approx.dfdu.isApprox(B, 1e-10));
  Eigen::Vector2d expectedF = A * x + B * u;
  EXPECT_TRUE(approx.f.isApprox(expectedF, 1e-10));
}

// 验证受控流映射使用控制器插值得到的输入。
TEST(DynamicsTest, ControlledFlowMap_UsesInterpolatedControllerInput) {
  Eigen::Matrix2d A = Eigen::Matrix2d::Zero();
  Eigen::Matrix2d B = Eigen::Matrix2d::Identity();
  LinearSystemDynamics<double, 2, 2> dyn(A, B);

  LinearController<double, 2, 2, 2> controller;
  controller.timeStamp_[0] = 0.0;
  controller.timeStamp_[1] = 1.0;
  controller.biasArray_[0] << 1.0, 2.0;
  controller.biasArray_[1] << 3.0, 4.0;
  controller.gainArray_[0].setZero();
  controller.gainArray_[1].setZero();
  dyn.setController(&controller);

  Eigen::Vector2d x;
  x << 10.0, -10.0;

  Eigen::Vector2d expected;
  expected << 2.0, 3.0;

  EXPECT_TRUE(dyn.computeFlowMap(0.5, x).isApprox(expected, 1e-10));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
