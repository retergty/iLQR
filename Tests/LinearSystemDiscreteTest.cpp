#include <iostream>
#include "iLQR.hpp"
#include <gtest/gtest.h>
#include "LinearSystemDynamics.hpp"
#include "OperatingPoints.hpp"
#include "math.h"

TEST(LinearSystemTest, DiscreteEK2Test)
{
    Eigen::Matrix<double, 3, 3> A = Eigen::Matrix<double, 3, 3>::Random();
    Eigen::Matrix<double, 3, 2> B = Eigen::Matrix<double, 3, 2>::Random();

    LinearSystemDynamics<double, 3, 2> lin_sys(A, B);

    TimeTriggeredRollout<double, 3, 2> rollout(&lin_sys, 1);

    const double initTime = 0;
    const double finalTime = 10;
    Eigen::Vector<double, 3> initState = Eigen::Vector<double, 3>::Random();

    LinearController<double, 3, 2, 11> controller;

    for (int i = 0; i < 11; ++i)
    {
        controller.timeStamp_[i] = i;
        controller.biasArray_[i] = Eigen::Vector2d::Random();
        controller.gainArray_[i] = Eigen::Matrix<double, 2, 3>::Random();
    }

    std::array<double, 11> rolloutTimeTrajectory;
    std::array<Eigen::Vector3d, 11> rolloutStateTrajectory;
    std::array<Eigen::Vector2d, 11> rolloutInputTrajectory;

    RolloutTrajectoryPointer<double, 3, 2> rolloutTrajectoryPointer(rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(), rolloutInputTrajectory.data(), 11);
    rollout.run(initTime, initState, finalTime, &controller, rolloutTrajectoryPointer);

    EK2DynamicsDiscretizer<double, 3, 2> ek2Discretizer;

    std::array<double, 11> testTimeTrajectory = rolloutTimeTrajectory;
    std::array<Eigen::Vector3d, 11> testStateTrajectory;

    testStateTrajectory[0] = initState;

    for (int i = 0; i < 10; ++i)
    {
        testStateTrajectory[i + 1] = ek2Discretizer.discretize(lin_sys, testTimeTrajectory[i], testStateTrajectory[i], rolloutInputTrajectory[i], 1.0);
    }

    for (int i = 0; i < 11; ++i)
    {
        EXPECT_FLOAT_EQ(testTimeTrajectory[i], rolloutTimeTrajectory[i]);
        EXPECT_TRUE(testStateTrajectory[i].isApprox(rolloutStateTrajectory[i], 1)) << "i = " << i;
    }
}

TEST(LinearSystemTest, DiscreteEK4Test)
{
    Eigen::Matrix<double, 3, 3> A = Eigen::Matrix<double, 3, 3>::Random();
    Eigen::Matrix<double, 3, 2> B = Eigen::Matrix<double, 3, 2>::Random();

    LinearSystemDynamics<double, 3, 2> lin_sys(A, B);

    TimeTriggeredRollout<double, 3, 2> rollout(&lin_sys, 1);

    const double initTime = 0;
    const double finalTime = 10;
    Eigen::Vector<double, 3> initState = Eigen::Vector<double, 3>::Random();

    LinearController<double, 3, 2, 11> controller;

    for (int i = 0; i < 11; ++i)
    {
        controller.timeStamp_[i] = i;
        controller.biasArray_[i] = Eigen::Vector2d::Zero();
        controller.gainArray_[i] = Eigen::Matrix<double, 2, 3>::Zero();
    }

    std::array<double, 11> rolloutTimeTrajectory;
    std::array<Eigen::Vector3d, 11> rolloutStateTrajectory;
    std::array<Eigen::Vector2d, 11> rolloutInputTrajectory;

    RolloutTrajectoryPointer<double, 3, 2> rolloutTrajectoryPointer(rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(), rolloutInputTrajectory.data(), 11);
    rollout.run(initTime, initState, finalTime, &controller, rolloutTrajectoryPointer);

    EK4DynamicsDiscretizer<double, 3, 2> ek4Discretizer;

    std::array<double, 11> testTimeTrajectory = rolloutTimeTrajectory;
    std::array<Eigen::Vector3d, 11> testStateTrajectory;

    testStateTrajectory[0] = initState;

    for (int i = 0; i < 10; ++i)
    {
        testStateTrajectory[i + 1] = ek4Discretizer.discretize(lin_sys, testTimeTrajectory[i], testStateTrajectory[i], rolloutInputTrajectory[i], 1.0);
    }

    for (int i = 0; i < 11; ++i)
    {
        EXPECT_FLOAT_EQ(testTimeTrajectory[i], rolloutTimeTrajectory[i]);
        EXPECT_TRUE(testStateTrajectory[i].isApprox(rolloutStateTrajectory[i], 0.1)) << "i = " << i;
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}