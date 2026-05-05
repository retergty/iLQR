#include <cmath>
#include <gtest/gtest.h>

#include "ModifiedRelaxedBarrierPenalty.hpp"
#include "Penalty.hpp"
#include "QuadraticPenalty.hpp"
#include "SlacknessSquaredHingePenalty.hpp"

// 验证二次罚函数的取值、导数和乘子更新公式。
TEST(PenaltyTest, QuadraticPenaltyValueDerivativesAndMultiplierUpdate)
{
    const QuadraticPenalty<double> penalty({4.0, 0.25});
    constexpr double time = 0.0;
    constexpr double lagrangian = 1.5;
    constexpr double constraint = -2.0;

    EXPECT_DOUBLE_EQ(penalty.getValue(time, lagrangian, constraint), 11.0);
    EXPECT_DOUBLE_EQ(penalty.getDerivative(time, lagrangian, constraint), -9.5);
    EXPECT_DOUBLE_EQ(penalty.getSecondDerivative(time, lagrangian, constraint), 4.0);
    EXPECT_DOUBLE_EQ(penalty.updateMultiplier(time, lagrangian, constraint), 3.5);
    EXPECT_DOUBLE_EQ(penalty.initializeMultiplier(), 0.0);
}

// 验证松弛平方 hinge 罚函数两个分支及其乘子更新。
TEST(PenaltyTest, SlacknessSquaredHingePenaltyBranchesAndMultiplierUpdate)
{
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
TEST(PenaltyTest, ModifiedRelaxedBarrierPenaltyMatchesLogAndRelaxedBranches)
{
    const ModifiedRelaxedBarrierPenalty<double> penalty({10.0, 0.1, 0.5});
    constexpr double time = 0.0;
    constexpr double lagrangian = 2.0;

    constexpr double logBranchConstraint = 0.05;
    EXPECT_NEAR(penalty.getValue(time, lagrangian, logBranchConstraint), -0.4 * std::log(1.25), 1e-12);
    EXPECT_NEAR(penalty.getDerivative(time, lagrangian, logBranchConstraint), -1.6, 1e-12);
    EXPECT_NEAR(penalty.getSecondDerivative(time, lagrangian, logBranchConstraint), 6.4, 1e-12);
    EXPECT_NEAR(penalty.updateMultiplier(time, lagrangian, logBranchConstraint), 1.6, 1e-12);

    constexpr double relaxedBranchConstraint = 0.0;
    constexpr double c2 = 1.0 / (1.1 * 1.1);
    constexpr double c1 = -1.0 / 1.1;
    const double c0 = -std::log(1.1);
    constexpr double vDelta = -0.1;
    const double expectedValue = 0.4 * (0.5 * c2 * vDelta * vDelta + c1 * vDelta + c0);
    const double expectedDerivative = 0.4 * (c2 * vDelta + c1) * 5.0;
    const double expectedSecondDerivative = 0.4 * c2 * 25.0;
    const double expectedMultiplier = 0.5 * 0.4 * (-c2 * vDelta - c1) * 5.0;

    EXPECT_NEAR(penalty.getValue(time, lagrangian, relaxedBranchConstraint), expectedValue, 1e-12);
    EXPECT_NEAR(penalty.getDerivative(time, lagrangian, relaxedBranchConstraint), expectedDerivative, 1e-12);
    EXPECT_NEAR(penalty.getSecondDerivative(time, lagrangian, relaxedBranchConstraint), expectedSecondDerivative, 1e-12);
    EXPECT_NEAR(penalty.updateMultiplier(time, lagrangian, relaxedBranchConstraint), expectedMultiplier, 1e-12);
    EXPECT_DOUBLE_EQ(penalty.initializeMultiplier(), 1.0);
}

// 验证罚函数包装器对纯状态线性约束正确应用链式法则。
TEST(PenaltyTest, PenaltyWrapperLinearStateConstraintUsesChainRule)
{
    QuadraticPenalty<double> quadraticPenalty({2.0, 0.5});
    Penalty<double, 2, 0> penalty(&quadraticPenalty);

    ScalarFunctionLinearApproximation<double, 2, 0> constraint;
    constraint.f = 2.0;
    constraint.dfdx << 1.0, -3.0;

    const auto approximation = penalty.getQuadraticApproximation(0.0, constraint, 0.5);

    Eigen::Vector2d expectedDfdx;
    expectedDfdx << 3.5, -10.5;
    Eigen::Matrix2d expectedDfdxx;
    expectedDfdxx << 2.0, -6.0,
                    -6.0, 18.0;

    EXPECT_DOUBLE_EQ(approximation.f, 3.0);
    EXPECT_TRUE(approximation.dfdx.isApprox(expectedDfdx, 1e-12));
    EXPECT_TRUE(approximation.dfdxx.isApprox(expectedDfdxx, 1e-12));
}

// 验证罚函数包装器对状态输入线性约束正确应用链式法则。
TEST(PenaltyTest, PenaltyWrapperLinearStateInputConstraintUsesChainRule)
{
    QuadraticPenalty<double> quadraticPenalty({4.0, 0.5});
    Penalty<double, 2, 1> penalty(&quadraticPenalty);

    ScalarFunctionLinearApproximation<double, 2, 1> constraint;
    constraint.f = 1.0;
    constraint.dfdx << 2.0, -1.0;
    constraint.dfdu << 0.5;

    const auto approximation = penalty.getQuadraticApproximation(0.0, constraint, 1.5);

    Eigen::Vector2d expectedDfdx;
    expectedDfdx << 5.0, -2.5;
    Eigen::Matrix<double, 1, 1> expectedDfdu;
    expectedDfdu << 1.25;
    Eigen::Matrix2d expectedDfdxx;
    expectedDfdxx << 16.0, -8.0,
                    -8.0, 4.0;
    Eigen::Matrix<double, 1, 2> expectedDfdux;
    expectedDfdux << 4.0, -2.0;
    Eigen::Matrix<double, 1, 1> expectedDfduu;
    expectedDfduu << 1.0;

    EXPECT_DOUBLE_EQ(approximation.f, 0.5);
    EXPECT_TRUE(approximation.dfdx.isApprox(expectedDfdx, 1e-12));
    EXPECT_TRUE(approximation.dfdu.isApprox(expectedDfdu, 1e-12));
    EXPECT_TRUE(approximation.dfdxx.isApprox(expectedDfdxx, 1e-12));
    EXPECT_TRUE(approximation.dfdux.isApprox(expectedDfdux, 1e-12));
    EXPECT_TRUE(approximation.dfduu.isApprox(expectedDfduu, 1e-12));
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
