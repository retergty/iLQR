/**
 * @file CostTest.cpp
 * @brief 代价模块测试：QuadraticStateCost、QuadraticStateInputCost 的值与二次近似。
 */
#include <gtest/gtest.h>
#include "QuadraticStateCost.hpp"
#include "Types.hpp"
#include <array>
#include <cmath>

// 验证状态代价在参考状态处取值为零。
TEST(CostTest, QuadraticStateCost_ValueAtReference)
{
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    QuadraticStateCost<double, 2, 3> cost(Q);

    std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> stateTraj;
    stateTraj[0] << 1.0, 0.0;
    stateTraj[1] << 2.0, 1.0;
    stateTraj[2] << 3.0, 2.0;

    double val = cost.getValue(1, stateTraj[1], timeTraj, stateTraj);
    EXPECT_DOUBLE_EQ(val, 0.0);
}

// 验证状态代价会惩罚相对参考轨迹的偏差。
TEST(CostTest, QuadraticStateCost_ValueDeviation)
{
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    QuadraticStateCost<double, 2, 3> cost(Q);

    std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> stateTraj;
    stateTraj[0] << 0.0, 0.0;
    stateTraj[1] << 0.0, 0.0;
    stateTraj[2] << 0.0, 0.0;

    Eigen::Vector2d x;
    x << 2.0, 0.0;
    double val = cost.getValue(1, x, timeTraj, stateTraj);
    EXPECT_DOUBLE_EQ(val, 0.5 * 4.0); // 0.5 * (2^2 + 0^2) = 2
}

// 验证状态代价的二次近似与解析形式一致。
TEST(CostTest, QuadraticStateCost_QuadraticApproximation)
{
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    QuadraticStateCost<double, 2, 3> cost(Q);

    std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> stateTraj;
    stateTraj[0].setZero();
    stateTraj[1].setZero();
    stateTraj[2].setZero();

    Eigen::Vector2d x;
    x << 1.0, 2.0;
    auto approx = cost.getQuadraticApproximation(1, x, timeTraj, stateTraj);

    EXPECT_TRUE(approx.dfdxx.isApprox(Q, 1e-10));
    EXPECT_TRUE(approx.dfdx.isApprox(x, 1e-10)); // Q*x = I*x = x
    EXPECT_DOUBLE_EQ(approx.f, 0.5 * (1.0 + 4.0)); // 0.5 * x'x = 2.5
}

// 验证状态代价在时间上使用插值后的参考状态。
TEST(CostTest, QuadraticStateCost_ValueUsesInterpolatedReference)
{
    Eigen::Matrix2d Q;
    Q << 2.0, 0.0,
         0.0, 4.0;
    QuadraticStateCost<double, 2, 3> cost(Q);

    std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> stateTraj;
    stateTraj[0] << 0.0, 0.0;
    stateTraj[1] << 2.0, 4.0;
    stateTraj[2] << 4.0, 8.0;

    Eigen::Vector2d x;
    x << 3.0, 1.0;

    const double val = cost.getValue(0.5, x, timeTraj, stateTraj);

    EXPECT_DOUBLE_EQ(val, 6.0);
}

// 验证只有 Q 和 R 项时状态输入代价的取值。
TEST(CostTest, QuadraticStateInputCost_ValueQROnly)
{
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d R = Eigen::Matrix2d::Identity();
    QuadraticStateInputCost<double, 2, 2, 3> cost(Q, R);

    std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> stateTraj, inputTraj;
    stateTraj[0].setZero();
    stateTraj[1].setZero();
    stateTraj[2].setZero();
    inputTraj[0].setZero();
    inputTraj[1].setZero();
    inputTraj[2].setZero();

    Eigen::Vector2d x, u;
    x << 1.0, 0.0;
    u << 0.0, 1.0;
    double val = cost.getValue(1, x, u, timeTraj, stateTraj, inputTraj);
    EXPECT_DOUBLE_EQ(val, 0.5 * 1.0 + 0.5 * 1.0); // 0.5*x'Q*x + 0.5*u'R*u = 0.5 + 0.5 = 1 (with P=0)
}

// 验证只有 Q 和 R 项时状态输入代价的二次近似。
TEST(CostTest, QuadraticStateInputCost_QuadraticApproximationQROnly)
{
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d R = 2.0 * Eigen::Matrix2d::Identity();
    QuadraticStateInputCost<double, 2, 2, 3> cost(Q, R);

    std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> stateTraj, inputTraj;
    stateTraj[0].setZero();
    stateTraj[1].setZero();
    stateTraj[2].setZero();
    inputTraj[0].setZero();
    inputTraj[1].setZero();
    inputTraj[2].setZero();

    Eigen::Vector2d x, u;
    x << 1.0, 1.0;
    u << 1.0, 0.0;
    auto approx = cost.getQuadraticApproximation(1, x, u, timeTraj, stateTraj, inputTraj);

    EXPECT_TRUE(approx.dfdxx.isApprox(Q, 1e-10));
    EXPECT_TRUE(approx.dfduu.isApprox(R, 1e-10));
    EXPECT_TRUE(approx.dfdx.isApprox(x, 1e-10));
    EXPECT_TRUE(approx.dfdu.isApprox(R * u, 1e-10));
    EXPECT_NEAR(approx.f, 0.5 * 2.0 + 0.5 * 2.0, 1e-10); // 0.5*x'x + 0.5*u'R*u = 1 + 1 = 2
}

// 验证状态输入交叉项会正确贡献到节点处的标量代价。
TEST(CostTest, QuadraticStateInputCost_ValueIncludesCrossTermAtIndex)
{
    Eigen::Matrix2d Q;
    Q << 2.0, 0.0,
         0.0, 4.0;
    Eigen::Matrix2d R;
    R << 3.0, 0.0,
         0.0, 5.0;
    Eigen::Matrix2d P;
    P << 0.5, -1.0,
         2.0, 0.25;
    QuadraticStateInputCost<double, 2, 2, 3> cost(Q, R, P);

    std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> stateTraj, inputTraj;
    stateTraj[0].setZero();
    stateTraj[1] << 1.0, 2.0;
    stateTraj[2].setZero();
    inputTraj[0].setZero();
    inputTraj[1] << 1.0, 1.0;
    inputTraj[2].setZero();

    Eigen::Vector2d x, u;
    x << 3.0, -1.0;
    u << 4.0, -2.0;

    const Eigen::Vector2d dx = x - stateTraj[1];
    const Eigen::Vector2d du = u - inputTraj[1];
    const double expected = 0.5 * dx.dot(Q * dx) + 0.5 * du.dot(R * du) + du.dot(P * dx);

    EXPECT_DOUBLE_EQ(cost.getValue(1, x, u, timeTraj, stateTraj, inputTraj), expected);
}

// 验证状态输入交叉项会正确出现在二次近似中。
TEST(CostTest, QuadraticStateInputCost_QuadraticApproximationIncludesCrossTerm)
{
    Eigen::Matrix2d Q;
    Q << 2.0, 0.0,
         0.0, 4.0;
    Eigen::Matrix2d R;
    R << 3.0, 0.0,
         0.0, 5.0;
    Eigen::Matrix2d P;
    P << 0.5, -1.0,
         2.0, 0.25;
    QuadraticStateInputCost<double, 2, 2, 3> cost(Q, R, P);

    std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> stateTraj, inputTraj;
    stateTraj[0] << 0.0, 0.0;
    stateTraj[1] << 2.0, 2.0;
    stateTraj[2] << 4.0, 4.0;
    inputTraj[0] << 0.0, 0.0;
    inputTraj[1] << 1.0, -1.0;
    inputTraj[2] << 2.0, -2.0;

    Eigen::Vector2d x, u;
    x << 2.0, 5.0;
    u << 4.0, -3.0;

    const Eigen::Vector2d referenceX = 0.5 * stateTraj[0] + 0.5 * stateTraj[1];
    const Eigen::Vector2d referenceU = 0.5 * inputTraj[0] + 0.5 * inputTraj[1];
    const Eigen::Vector2d dx = x - referenceX;
    const Eigen::Vector2d du = u - referenceU;

    const auto approx = cost.getQuadraticApproximation(0.5, x, u, timeTraj, stateTraj, inputTraj);

    EXPECT_TRUE(approx.dfdxx.isApprox(Q, 1e-10));
    EXPECT_TRUE(approx.dfduu.isApprox(R, 1e-10));
    EXPECT_TRUE(approx.dfdux.isApprox(P, 1e-10));
    EXPECT_TRUE(approx.dfdx.isApprox(Q * dx + P.transpose() * du, 1e-10));
    EXPECT_TRUE(approx.dfdu.isApprox(R * du + P * dx, 1e-10));
    EXPECT_NEAR(approx.f, 0.5 * dx.dot(Q * dx) + 0.5 * du.dot(R * du) + du.dot(P * dx), 1e-10);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
