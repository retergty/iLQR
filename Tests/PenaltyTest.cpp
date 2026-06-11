#include <gtest/gtest.h>

#include <cmath>

#include "Penalties/ModifiedRelaxedBarrierPenalty.hpp"
#include "Penalties/QuadraticPenalty.hpp"
#include "Penalties/SlacknessSquaredHingePenalty.hpp"

// 验证二次罚函数的取值、导数和乘子更新公式。
TEST(PenaltyTest, QuadraticPenaltyValueDerivativesAndMultiplierUpdate) {
  const QuadraticPenalty<double> penalty({4.0, 0.25});
  constexpr double time = 0.0;
  constexpr double lagrangian = 1.5;
  constexpr double constraint = -2.0;

  EXPECT_DOUBLE_EQ(penalty.getValue(time, lagrangian, constraint), 11.0);
  EXPECT_DOUBLE_EQ(penalty.getDerivative(time, lagrangian, constraint), -9.5);
  EXPECT_DOUBLE_EQ(penalty.getSecondDerivative(time, lagrangian, constraint),
                   4.0);
  EXPECT_DOUBLE_EQ(penalty.updateMultiplier(time, lagrangian, constraint), 3.5);
  EXPECT_DOUBLE_EQ(penalty.initializeMultiplier(), 0.0);
}

// 验证松弛平方 hinge 罚函数两个分支及其乘子更新。
TEST(PenaltyTest, SlacknessSquaredHingePenaltyBranchesAndMultiplierUpdate) {
  const SlacknessSquaredHingePenalty<double> penalty({10.0, 0.5});
  constexpr double time = 0.0;
  constexpr double lagrangian = 2.0;

  EXPECT_DOUBLE_EQ(penalty.getValue(time, lagrangian, 0.1), -0.15);
  EXPECT_DOUBLE_EQ(penalty.getDerivative(time, lagrangian, 0.1), -1.0);
  EXPECT_DOUBLE_EQ(penalty.getSecondDerivative(time, lagrangian, 0.1), 10.0);
  EXPECT_DOUBLE_EQ(penalty.updateMultiplier(time, lagrangian, 0.1), 1.5);

  EXPECT_DOUBLE_EQ(penalty.getValue(time, lagrangian, 1.0), -0.2);
  EXPECT_DOUBLE_EQ(penalty.getDerivative(time, lagrangian, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(penalty.getSecondDerivative(time, lagrangian, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(penalty.updateMultiplier(time, lagrangian, 1.0), 1.0);
  EXPECT_DOUBLE_EQ(penalty.initializeMultiplier(), 0.0);
}

// 验证改进松弛 barrier 罚函数与对数分支和松弛分支公式一致。
TEST(PenaltyTest, ModifiedRelaxedBarrierPenaltyMatchesLogAndRelaxedBranches) {
  const ModifiedRelaxedBarrierPenalty<double> penalty({10.0, 0.1, 0.5});
  constexpr double time = 0.0;
  constexpr double lagrangian = 2.0;

  constexpr double logBranchConstraint = 0.05;
  EXPECT_NEAR(penalty.getValue(time, lagrangian, logBranchConstraint),
              -0.4 * std::log(1.25), 1e-12);
  EXPECT_NEAR(penalty.getDerivative(time, lagrangian, logBranchConstraint),
              -1.6, 1e-12);
  EXPECT_NEAR(
      penalty.getSecondDerivative(time, lagrangian, logBranchConstraint), 6.4,
      1e-12);
  EXPECT_NEAR(penalty.updateMultiplier(time, lagrangian, logBranchConstraint),
              1.6, 1e-12);

  constexpr double relaxedBranchConstraint = 0.0;
  constexpr double c2 = 1.0 / (1.1 * 1.1);
  constexpr double c1 = -1.0 / 1.1;
  const double c0 = -std::log(1.1);
  constexpr double vDelta = -0.1;
  const double expectedValue =
      0.4 * (0.5 * c2 * vDelta * vDelta + c1 * vDelta + c0);
  const double expectedDerivative = 0.4 * (c2 * vDelta + c1) * 5.0;
  const double expectedSecondDerivative = 0.4 * c2 * 25.0;
  const double expectedMultiplier = 0.5 * 0.4 * (-c2 * vDelta - c1) * 5.0;

  EXPECT_NEAR(penalty.getValue(time, lagrangian, relaxedBranchConstraint),
              expectedValue, 1e-12);
  EXPECT_NEAR(penalty.getDerivative(time, lagrangian, relaxedBranchConstraint),
              expectedDerivative, 1e-12);
  EXPECT_NEAR(
      penalty.getSecondDerivative(time, lagrangian, relaxedBranchConstraint),
      expectedSecondDerivative, 1e-12);
  EXPECT_NEAR(
      penalty.updateMultiplier(time, lagrangian, relaxedBranchConstraint),
      expectedMultiplier, 1e-12);
  EXPECT_DOUBLE_EQ(penalty.initializeMultiplier(), 1.0);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
