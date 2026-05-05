#include <iostream>
#include "iLQR.hpp"
#include <gtest/gtest.h>
#include "LinearSystemDynamics.hpp"
#include "OperatingPoints.hpp"
#include "math.h"

TEST(LinearSystemTest, DiscreteEK2Test)
{
    Eigen::Matrix2d A;
    A << 0.0, 1.0,
         0.0, 0.0;

    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;

    LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
    EK2DynamicsDiscretizer<double, 2, 1> ek2Discretizer;

    Eigen::Vector2d x;
    x << 1.0, 2.0;

    Eigen::Vector<double, 1> u;
    u << 3.0;

    const double t = 0.0;
    const double dt = 0.5;
    const Eigen::Vector2d k1 = A * x + B * u;
    const Eigen::Vector2d k2 = A * (x + dt * k1) + B * u;
    const Eigen::Vector2d expected = x + 0.5 * dt * (k1 + k2);

    const Eigen::Vector2d actual = ek2Discretizer.discretize(lin_sys, t, x, u, dt);

    EXPECT_TRUE(actual.isApprox(expected, 1e-12));
}

TEST(LinearSystemTest, DiscreteEulerTest)
{
    Eigen::Matrix2d A;
    A << 1.0, 2.0,
         0.0, -1.0;

    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 2.0;

    LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
    EulerDynamicsDiscretizer<double, 2, 1> eulerDiscretizer;

    Eigen::Vector2d x;
    x << 1.0, -2.0;

    Eigen::Vector<double, 1> u;
    u << 3.0;

    const double dt = 0.25;
    const Eigen::Vector2d expected = x + dt * (A * x + B * u);
    const Eigen::Vector2d actual = eulerDiscretizer.discretize(lin_sys, 1.0, x, u, dt);

    EXPECT_TRUE(actual.isApprox(expected, 1e-12));
}

TEST(LinearSystemTest, EulerSensitivityMatchesLinearSystemAnalyticResult)
{
    Eigen::Matrix2d A;
    A << 1.0, 2.0,
         0.0, -1.0;

    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 2.0;

    LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
    EulerDynamicsDiscretizer<double, 2, 1> eulerDiscretizer;

    Eigen::Vector2d x;
    x << 1.0, -2.0;

    Eigen::Vector<double, 1> u;
    u << 3.0;

    const double dt = 0.25;
    const auto approx = eulerDiscretizer.sensitivityDiscretize(lin_sys, 1.0, x, u, dt);

    EXPECT_TRUE(approx.dfdx.isApprox(Eigen::Matrix2d::Identity() + dt * A, 1e-12));
    EXPECT_TRUE(approx.dfdu.isApprox(dt * B, 1e-12));
    EXPECT_TRUE(approx.f.isApprox(x + dt * (A * x + B * u), 1e-12));
}

TEST(LinearSystemTest, EK2SensitivityMatchesLinearSystemAnalyticResult)
{
    Eigen::Matrix2d A;
    A << 0.0, 1.0,
         -2.0, -3.0;

    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;

    LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
    EK2DynamicsDiscretizer<double, 2, 1> ek2Discretizer;

    Eigen::Vector2d x;
    x << 1.0, -1.0;

    Eigen::Vector<double, 1> u;
    u << 2.0;

    const double dt = 0.1;
    const auto approx = ek2Discretizer.sensitivityDiscretize(lin_sys, 0.0, x, u, dt);
    const Eigen::Matrix2d expectedDfdx = Eigen::Matrix2d::Identity() + dt * A + 0.5 * dt * dt * A * A;
    const Eigen::Matrix<double, 2, 1> expectedDfdu = dt * B + 0.5 * dt * dt * A * B;
    const Eigen::Vector2d expectedF = ek2Discretizer.discretize(lin_sys, 0.0, x, u, dt);

    EXPECT_TRUE(approx.dfdx.isApprox(expectedDfdx, 1e-12));
    EXPECT_TRUE(approx.dfdu.isApprox(expectedDfdu, 1e-12));
    EXPECT_TRUE(approx.f.isApprox(expectedF, 1e-12));
}

TEST(LinearSystemTest, DiscreteEK4Test)
{
    Eigen::Matrix2d A;
    A << 0.0, 1.0,
         -2.0, -3.0;

    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;

    LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
    EK4DynamicsDiscretizer<double, 2, 1> ek4Discretizer;

    Eigen::Vector2d x;
    x << 1.0, -1.0;

    Eigen::Vector<double, 1> u;
    u << 2.0;

    const double t = 0.0;
    const double dt = 0.1;
    const Eigen::Vector2d k1 = A * x + B * u;
    const Eigen::Vector2d k2 = A * (x + 0.5 * dt * k1) + B * u;
    const Eigen::Vector2d k3 = A * (x + 0.5 * dt * k2) + B * u;
    const Eigen::Vector2d k4 = A * (x + dt * k3) + B * u;
    const Eigen::Vector2d expected = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

    const Eigen::Vector2d actual = ek4Discretizer.discretize(lin_sys, t, x, u, dt);

    EXPECT_TRUE(actual.isApprox(expected, 1e-12));
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}