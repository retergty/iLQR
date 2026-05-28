#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "Controller/LinearController.hpp"
#include "Dynamics/LinearSystemDynamics.hpp"
#include "Rollout/TimeTriggeredRollout.hpp"

// 验证零输入 rollout 与 xdot = x 的解析解一致。
TEST(LinearSystemRolloutTest, ZeroInputMatchesExponentialSolution) {
  Matrix<double, 3, 3> A = Matrix<double, 3, 3>::Identity();
  Matrix<double, 3, 2> B{{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
  LinearSystemDynamics<double, 3, 2> lin_sys(A, B);

  TimeTriggeredRollout<double, 3, 2> rollout(&lin_sys, 0.1);

  const double initTime = 0;
  const double finalTime = 1;
  LinearController<double, 3, 2, 11> controller;

  for (int i = 0; i < 11; ++i) {
    controller.timeStamp_[i] = 0.1 * i;
    controller.biasArray_[i] = Vector<double, 2>::Zero();
    controller.gainArray_[i] = Matrix<double, 2, 3>::Zero();
  }

  std::array<double, 11> rolloutTimeTrajectory;
  std::array<Vector<double, 3>, 11> rolloutStateTrajectory;
  std::array<Vector<double, 2>, 11> rolloutInputTrajectory;

  RolloutTrajectoryPointer<double, 3, 2> rolloutTrajectoryPointer(
      rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(),
      rolloutInputTrajectory.data(), 11);
  const int count = rollout.run(initTime, Vector<double, 3>::Ones(), finalTime,
                                &controller, rolloutTrajectoryPointer);

  ASSERT_EQ(count, 11);
  for (int i = 0; i < count; ++i) {
    const double t = 0.1 * i;
    const Vector<double, 3> expectedState =
        Vector<double, 3>::Constant(std::exp(t));

    EXPECT_NEAR(rolloutTimeTrajectory[i], t, 1e-12);
    EXPECT_TRUE(rolloutStateTrajectory[i].isApprox(expectedState, 1e-6))
        << "i = " << i;
    EXPECT_TRUE(
        rolloutInputTrajectory[i].isApprox(Vector<double, 2>::Zero(), 1e-12))
        << "i = " << i;
  }
}

// 验证带反馈控制器的 rollout 与闭环解析解一致。
TEST(LinearSystemRolloutTest, FeedbackControllerMatchesClosedLoopSolution) {
  Matrix<double, 3, 3> A = Matrix<double, 3, 3>::Identity();
  Matrix<double, 3, 3> B = Matrix<double, 3, 3>::Identity();

  LinearSystemDynamics<double, 3, 3> lin_sys(A, B);

  TimeTriggeredRollout<double, 3, 3> rollout(&lin_sys, 0.1);

  const double initTime = 0;
  const double finalTime = 1;
  LinearController<double, 3, 3, 11> controller;

  for (int i = 0; i < 11; ++i) {
    controller.timeStamp_[i] = 0.1 * i;
    controller.biasArray_[i] = Vector<double, 3>::Zero();
    controller.gainArray_[i] = -2.0 * Matrix<double, 3, 3>::Identity();
  }

  std::array<double, 11> rolloutTimeTrajectory;
  std::array<Vector<double, 3>, 11> rolloutStateTrajectory;
  std::array<Vector<double, 3>, 11> rolloutInputTrajectory;

  RolloutTrajectoryPointer<double, 3, 3> rolloutTrajectoryPointer(
      rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(),
      rolloutInputTrajectory.data(), 11);
  const int count =
      rollout.run(initTime, 10.0 * Vector<double, 3>::Ones(), finalTime,
                  &controller, rolloutTrajectoryPointer);

  ASSERT_EQ(count, 11);
  for (int i = 0; i < count; ++i) {
    const double t = 0.1 * i;
    const Vector<double, 3> expectedState =
        Vector<double, 3>::Constant(10 * std::exp(-t));
    const Vector<double, 3> expectedInput = -2.0 * expectedState;

    EXPECT_NEAR(rolloutTimeTrajectory[i], t, 1e-12);
    EXPECT_TRUE(rolloutStateTrajectory[i].isApprox(expectedState, 1e-6))
        << "i = " << i;
    EXPECT_TRUE(rolloutInputTrajectory[i].isApprox(expectedInput, 1e-6))
        << "i = " << i;
  }
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}