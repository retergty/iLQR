#include <gtest/gtest.h>

#include <array>

#include "Dynamics/LinearSystemDynamics.hpp"
#include "iLQR/iLQR.hpp"

// 验证 incrementController 保留时间戳和增益，只更新前馈偏置。
TEST(LinearSystemILQRTest, IncrementControllerUpdatesOnlyFeedforwardBias) {
  using Descriptor =
      iLQRDescriptor<double, TranscriptionConfig<Dimensions<2, 1>, Horizon<2>>>;
  using Solver = iLQR<Descriptor>;
  Solver::LinearController_t unoptimizedController;
  Solver::LinearController_t controller;

  for (int i = 0; i < 3; ++i) {
    unoptimizedController.timeStamp_[i] = 0.5 * i;
    unoptimizedController.biasArray_[i] = {1.0 + i};
    unoptimizedController.deltaBiasArray_[i] = {4.0 - i};
    unoptimizedController.gainArray_[i] = {2.0 + i, -1.0 * i};

    controller.timeStamp_[i] = -1.0;
    controller.biasArray_[i] = {-100.0};
    controller.deltaBiasArray_[i] = {-200.0};
    controller.gainArray_[i] = Matrix<double, 1, 2>::Constant(-300.0);
  }

  Solver::incrementController(0.25, unoptimizedController, controller);

  for (int i = 0; i < 3; ++i) {
    const Vector<double, 1> expectedBias{(1.0 + i) + 0.25 * (4.0 - i)};

    EXPECT_DOUBLE_EQ(controller.timeStamp_[i],
                     unoptimizedController.timeStamp_[i]);
    EXPECT_TRUE(controller.gainArray_[i].isApprox(
        unoptimizedController.gainArray_[i], 1e-12))
        << "i = " << i;
    EXPECT_TRUE(controller.biasArray_[i].isApprox(expectedBias, 1e-12))
        << "i = " << i;
  }
}

// 验证控制器更新量的积分按梯形积分计算。
TEST(LinearSystemILQRTest, ControllerUpdateIntegralUsesTrapezoidalRule) {
  using Descriptor =
      iLQRDescriptor<double, TranscriptionConfig<Dimensions<2, 1>, Horizon<2>>>;
  using Solver = iLQR<Descriptor>;
  Solver::LinearController_t controller;

  controller.timeStamp_[0] = 0.0;
  controller.timeStamp_[1] = 1.0;
  controller.timeStamp_[2] = 3.0;
  controller.deltaBiasArray_[0] = {1.0};
  controller.deltaBiasArray_[1] = {3.0};
  controller.deltaBiasArray_[2] = {5.0};

  // 在 [0, 1, 3] 上对平方范数 [1, 9, 25] 积分。
  EXPECT_DOUBLE_EQ(Solver::computeControllerUpdateIS(controller), 39.0);
}

// 验证 rolloutTrajectory 会把时间、状态和输入写回 primal solution。
TEST(LinearSystemILQRTest, RolloutTrajectoryWritesPrimalSolution) {
  using Descriptor =
      iLQRDescriptor<double, TranscriptionConfig<Dimensions<2, 2>, Horizon<4>>>;
  using Solver = iLQR<Descriptor>;

  Matrix<double, 2, 2> A = Matrix<double, 2, 2>::Zero();
  Matrix<double, 2, 2> B = Matrix<double, 2, 2>::Identity();
  LinearSystemDynamics<double, 2, 2> dynamics(A, B);
  TimeTriggeredRollout<double, 2, 2> rollout(&dynamics, 0.5);
  Solver::PrimalSolution_t primalSolution;

  for (int i = 0; i < 5; ++i) {
    primalSolution.controller_.timeStamp_[i] = 0.5 * i;
    primalSolution.controller_.biasArray_[i] = {1.0, -2.0};
    primalSolution.controller_.gainArray_[i] = Matrix<double, 2, 2>::Zero();
  }

  const Vector<double, 2> initState{3.0, 4.0};
  const double averageTimeStep =
      Solver::rolloutTrajectory(rollout, 0.0, initState, 2.0, primalSolution);

  EXPECT_DOUBLE_EQ(averageTimeStep, 0.5);
  for (int i = 0; i < 5; ++i) {
    const double t = 0.5 * i;
    const Vector<double, 2> expectedState{3.0 + t, 4.0 - 2.0 * t};
    const Vector<double, 2> expectedInput{1.0, -2.0};

    EXPECT_DOUBLE_EQ(primalSolution.timeTrajectory_[i], t);
    EXPECT_TRUE(
        primalSolution.stateTrajectory_[i].isApprox(expectedState, 1e-10))
        << "i = " << i;
    EXPECT_TRUE(
        primalSolution.inputTrajectory_[i].isApprox(expectedInput, 1e-10))
        << "i = " << i;
  }
}

// 验证在 A = 0、B = I 时，常值开环输入的 rollout 与解析解一致。
TEST(TimeTriggeredRolloutTest, ConstantInputMatchesAnalyticSolution) {
  Matrix<double, 2, 2> A = Matrix<double, 2, 2>::Zero();
  Matrix<double, 2, 2> B = Matrix<double, 2, 2>::Identity();
  LinearSystemDynamics<double, 2, 2> lin_sys(A, B);

  TimeTriggeredRollout<double, 2, 2> rollout(&lin_sys, 0.5);
  LinearController<double, 2, 2, 5> controller;
  for (int i = 0; i < 5; ++i) {
    controller.timeStamp_[i] = 0.5 * i;
    controller.biasArray_[i] = {1.0, -2.0};
    controller.gainArray_[i] = Matrix<double, 2, 2>::Zero();
  }

  std::array<double, 5> rolloutTimeTrajectory;
  std::array<Vector<double, 2>, 5> rolloutStateTrajectory;
  std::array<Vector<double, 2>, 5> rolloutInputTrajectory;
  RolloutTrajectoryPointer<double, 2, 2> rolloutTrajectoryPointer(
      rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(),
      rolloutInputTrajectory.data(), 5);

  const Vector<double, 2> initState{3.0, 4.0};
  const int count =
      rollout.run(0.0, initState, 2.0, &controller, rolloutTrajectoryPointer);

  EXPECT_EQ(count, 5);
  for (int i = 0; i < count; ++i) {
    const double t = 0.5 * i;
    const Vector<double, 2> expectedState{3.0 + t, 4.0 - 2.0 * t};
    const Vector<double, 2> expectedInput{1.0, -2.0};

    EXPECT_DOUBLE_EQ(rolloutTimeTrajectory[i], t);
    EXPECT_TRUE(rolloutStateTrajectory[i].isApprox(expectedState, 1e-10))
        << "i = " << i;
    EXPECT_TRUE(rolloutInputTrajectory[i].isApprox(expectedInput, 1e-10))
        << "i = " << i;
  }
}

// 验证 rollout 的时间戳和状态演化正确处理非零初始时刻。
TEST(TimeTriggeredRolloutTest, RespectsNonZeroInitialTime) {
  Matrix<double, 2, 2> A = Matrix<double, 2, 2>::Zero();
  Matrix<double, 2, 2> B = Matrix<double, 2, 2>::Identity();
  LinearSystemDynamics<double, 2, 2> lin_sys(A, B);

  TimeTriggeredRollout<double, 2, 2> rollout(&lin_sys, 0.25);
  LinearController<double, 2, 2, 5> controller;
  for (int i = 0; i < 5; ++i) {
    controller.timeStamp_[i] = 1.0 + 0.25 * i;
    controller.biasArray_[i] = {-2.0, 1.0};
    controller.gainArray_[i] = Matrix<double, 2, 2>::Zero();
  }

  std::array<double, 5> rolloutTimeTrajectory;
  std::array<Vector<double, 2>, 5> rolloutStateTrajectory;
  std::array<Vector<double, 2>, 5> rolloutInputTrajectory;
  RolloutTrajectoryPointer<double, 2, 2> rolloutTrajectoryPointer(
      rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(),
      rolloutInputTrajectory.data(), 5);

  const Vector<double, 2> initState{10.0, -10.0};
  const int count =
      rollout.run(1.0, initState, 2.0, &controller, rolloutTrajectoryPointer);

  EXPECT_EQ(count, 5);
  for (int i = 0; i < count; ++i) {
    const double t = 1.0 + 0.25 * i;
    const double elapsed = t - 1.0;
    const Vector<double, 2> expectedState{10.0 - 2.0 * elapsed,
                                          -10.0 + elapsed};

    EXPECT_DOUBLE_EQ(rolloutTimeTrajectory[i], t);
    EXPECT_TRUE(rolloutStateTrajectory[i].isApprox(expectedState, 1e-10))
        << "i = " << i;
  }
}

// 验证关闭输入重建后，调用方提供的输入缓冲区不会被改写。
TEST(TimeTriggeredRolloutTest, CanSkipInputTrajectoryReconstruction) {
  Matrix<double, 1, 1> A{0.0};
  Matrix<double, 1, 1> B{1.0};
  LinearSystemDynamics<double, 1, 1> lin_sys(A, B);

  TimeTriggeredRollout<double, 1, 1> rollout(&lin_sys, 0.5);
  rollout.settings().reconstructInputTrajectory = false;

  LinearController<double, 1, 1, 3> controller;
  for (int i = 0; i < 3; ++i) {
    controller.timeStamp_[i] = 0.5 * i;
    controller.biasArray_[i] = {2.0};
    controller.gainArray_[i] = Matrix<double, 1, 1>::Zero();
  }

  std::array<double, 3> rolloutTimeTrajectory;
  std::array<Vector<double, 1>, 3> rolloutStateTrajectory;
  std::array<Vector<double, 1>, 3> rolloutInputTrajectory;
  for (auto& input : rolloutInputTrajectory) {
    input = {-123.0};
  }

  RolloutTrajectoryPointer<double, 1, 1> rolloutTrajectoryPointer(
      rolloutTimeTrajectory.data(), rolloutStateTrajectory.data(),
      rolloutInputTrajectory.data(), 3);

  const Vector<double, 1> initState{1.0};
  const int count =
      rollout.run(0.0, initState, 1.0, &controller, rolloutTrajectoryPointer);

  EXPECT_EQ(count, 3);
  for (int i = 0; i < count; ++i) {
    EXPECT_DOUBLE_EQ(rolloutTimeTrajectory[i], 0.5 * i);
    EXPECT_DOUBLE_EQ(rolloutStateTrajectory[i](0), 1.0 + i);
    EXPECT_DOUBLE_EQ(rolloutInputTrajectory[i](0), -123.0);
  }
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}