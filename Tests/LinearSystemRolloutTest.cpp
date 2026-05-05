#include <array>
#include <cmath>
#include <gtest/gtest.h>

#include "TimeTriggeredRollout.hpp"
#include "LinearController.hpp"
#include "LinearSystemDynamics.hpp"

// 验证零输入 rollout 与 xdot = x 的解析解一致。
TEST(LinearSystemRolloutTest, ZeroInputMatchesExponentialSolution)
{
    Eigen::Matrix<double, 3, 3> A = Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Matrix<double, 3, 2> B;
    B << 1, 0,
        0, 1,
        1, 1;
    LinearSystemDynamics<double, 3, 2> lin_sys(A, B);

    TimeTriggeredRollout<double, 3, 2> rollout(&lin_sys, 0.1);

    const double initTime = 0;
    const double finalTime = 1;
    LinearController<double, 3, 2, 11> controller;

    for (int i = 0; i < 11; ++i)
    {
        controller.timeStamp_[i] = 0.1 * i;
        controller.biasArray_[i] = Eigen::Vector2d::Zero();
        controller.gainArray_[i] = Eigen::Matrix<double, 2, 3>::Zero();
    }

    std::array<double, 11> rolloutTimeTrajectory;
    std::array<Eigen::Vector3d, 11> rolloutStateTrajectory;
    std::array<Eigen::Vector2d, 11> rolloutInputTrajectory;

    RolloutTrajectoryPointer<double, 3, 2> rolloutTrajectoryPointer(rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(), rolloutInputTrajectory.data(), 11);
    const int count = rollout.run(initTime, Eigen::Vector3d::Ones(), finalTime, &controller, rolloutTrajectoryPointer);

    ASSERT_EQ(count, 11);
    for (int i = 0; i < count; ++i)
    {
        const double t = 0.1 * i;
        Eigen::Vector3d expectedState;
        expectedState.setConstant(std::exp(t));

        EXPECT_NEAR(rolloutTimeTrajectory[i], t, 1e-12);
        EXPECT_TRUE(rolloutStateTrajectory[i].isApprox(expectedState, 1e-6)) << "i = " << i;
        EXPECT_TRUE(rolloutInputTrajectory[i].isApprox(Eigen::Vector2d::Zero(), 1e-12)) << "i = " << i;
    }
}

// 验证带反馈控制器的 rollout 与闭环解析解一致。
TEST(LinearSystemRolloutTest, FeedbackControllerMatchesClosedLoopSolution)
{
    Eigen::Matrix<double, 3, 3> A = Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Matrix<double, 3, 3> B = Eigen::Matrix<double, 3, 3>::Identity();

    LinearSystemDynamics<double, 3, 3> lin_sys(A, B);

    TimeTriggeredRollout<double, 3, 3> rollout(&lin_sys, 0.1);

    const double initTime = 0;
    const double finalTime = 1;
    LinearController<double, 3, 3, 11> controller;

    for (int i = 0; i < 11; ++i)
    {
        controller.timeStamp_[i] = 0.1 * i;
        controller.biasArray_[i] = Eigen::Vector3d::Zero();
        controller.gainArray_[i] = -2 * Eigen::Matrix<double, 3, 3>::Identity();
    }

    std::array<double, 11> rolloutTimeTrajectory;
    std::array<Eigen::Vector3d, 11> rolloutStateTrajectory;
    std::array<Eigen::Vector3d, 11> rolloutInputTrajectory;

    RolloutTrajectoryPointer<double, 3, 3> rolloutTrajectoryPointer(rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(), rolloutInputTrajectory.data(), 11);
    const int count = rollout.run(initTime, 10 * Eigen::Vector3d::Ones(), finalTime, &controller, rolloutTrajectoryPointer);

    ASSERT_EQ(count, 11);
    for (int i = 0; i < count; ++i)
    {
        const double t = 0.1 * i;
        Eigen::Vector3d expectedState;
        expectedState.setConstant(10 * std::exp(-t));
        const Eigen::Vector3d expectedInput = -2 * expectedState;

        EXPECT_NEAR(rolloutTimeTrajectory[i], t, 1e-12);
        EXPECT_TRUE(rolloutStateTrajectory[i].isApprox(expectedState, 1e-6)) << "i = " << i;
        EXPECT_TRUE(rolloutInputTrajectory[i].isApprox(expectedInput, 1e-6)) << "i = " << i;
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}