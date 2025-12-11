#include <iostream>
#include "iLQR.hpp"
#include <gtest/gtest.h>
#include "LinearSystemDynamics.hpp"
#include "OperatingPoints.hpp"
#include "math.h"
TEST(LinearSystemTest, Rollout)
{
  Eigen::Matrix<double, 3, 3> A = Eigen::Matrix<double, 3, 3>::Identity();
  Eigen::Matrix<double, 3, 2> B;
  B << 1, 0,
    0, 1,
    1, 1;
  LinearSystemDynamics<double, 3, 2> lin_sys(A, B);

  TimeTriggeredRollout<double, 3, 2> rollout(&lin_sys, 1);

  const double initTime = 0;
  const double finalTime = 10;
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
  rollout.run(initTime, Eigen::Vector3d::Ones(), finalTime, &controller, rolloutTrajectoryPointer);

  std::array<double, 11> testTimeTrajectory;
  std::array<Eigen::Vector3d, 11> testStateTrajectory;
  std::array<Eigen::Vector2d, 11> testInputTrajectory;

  for (int i = 0; i < 11; ++i)
  {
    testTimeTrajectory[i] = i;
    testStateTrajectory[i].setConstant(std::exp(i));
    testInputTrajectory[i].setZero();
  }

  for (int i = 0; i < 11; ++i)
  {
    EXPECT_FLOAT_EQ(testTimeTrajectory[i], rolloutTimeTrajectory[i]);
    EXPECT_TRUE(testStateTrajectory[i].isApprox(rolloutStateTrajectory[i], 1)) << "i = " << i;
    EXPECT_TRUE(testInputTrajectory[i].isApprox(rolloutInputTrajectory[i], 1)) << "i = " << i;
  }
}

TEST(LinearSystemTest, Rollout2)
{
  Eigen::Matrix<double, 3, 3> A = Eigen::Matrix<double, 3, 3>::Identity();
  Eigen::Matrix<double, 3, 3> B = Eigen::Matrix<double, 3, 3>::Identity();

  LinearSystemDynamics<double, 3, 3> lin_sys(A, B);

  TimeTriggeredRollout<double, 3, 3> rollout(&lin_sys, 1);

  const double initTime = 0;
  const double finalTime = 10;
  LinearController<double, 3, 3, 11> controller;

  for (int i = 0; i < 11; ++i)
  {
    controller.timeStamp_[i] = i;
    controller.biasArray_[i] = Eigen::Vector3d::Zero();
    controller.gainArray_[i] = -2 * Eigen::Matrix<double, 3, 3>::Identity();
  }

  std::array<double, 11> rolloutTimeTrajectory;
  std::array<Eigen::Vector3d, 11> rolloutStateTrajectory;
  std::array<Eigen::Vector3d, 11> rolloutInputTrajectory;

  RolloutTrajectoryPointer<double, 3, 3> rolloutTrajectoryPointer(rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(), rolloutInputTrajectory.data(), 11);
  rollout.run(initTime, 10 * Eigen::Vector3d::Ones(), finalTime, &controller, rolloutTrajectoryPointer);

  std::array<double, 11> testTimeTrajectory;
  std::array<Eigen::Vector3d, 11> testStateTrajectory;
  std::array<Eigen::Vector3d, 11> testInputTrajectory;

  for (int i = 0; i < 11; ++i)
  {
    testTimeTrajectory[i] = i;
    testStateTrajectory[i].setConstant(10 * std::exp(-i));
    testInputTrajectory[i].noalias() = -2 * Eigen::Matrix3d::Identity() * testStateTrajectory[i];
  }

  for (int i = 0; i < 11; ++i)
  {
    EXPECT_FLOAT_EQ(testTimeTrajectory[i], rolloutTimeTrajectory[i]);
    EXPECT_TRUE((testStateTrajectory[i] - rolloutStateTrajectory[i]).norm() < 1e-1) << "i = " << i;
    EXPECT_TRUE((testInputTrajectory[i] - rolloutInputTrajectory[i]).norm() < 1e-1) << "i = " << i;
  }
}



int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}