/**
 * @file iLQREndToEndTest.cpp
 * @brief iLQR 求解器端到端测试：构造问题并运行求解。
 */
#include <gtest/gtest.h>

#include "DefaultInitializer.hpp"
#include "LinearSystemDynamics.hpp"
#include "QuadraticStateCost.hpp"
#include "Types.hpp"
#include "iLQR.hpp"

// 验证带二次跟踪代价的问题可以完整运行 iLQR 求解流程。
TEST(iLQREndToEndTest, RunWithQuadraticTrackingCost) {
  Eigen::Matrix2d A = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d B = Eigen::Matrix2d::Identity();
  LinearSystemDynamics<double, 2, 2> dynamics(A, B);
  DefaultInitializer<double, 2, 2> initializer;
  DDPSettings<double> ddp_setting;
  using Descriptor =
      iLQRDescriptor<double, TranscriptionConfig<Dimensions<2, 2>, Horizon<5>>>;
  iLQR<Descriptor> solver(ddp_setting, &dynamics, &initializer);
  QuadraticStateInputCost<double, 2, 2, 6> runningCost(
      Eigen::Matrix2d::Identity(), Eigen::Matrix2d::Identity());
  QuadraticStateCost<double, 2, 6> finalCost(2.0 * Eigen::Matrix2d::Identity());
  solver.optimalControlProblem_.cost.add(runningCost);
  solver.optimalControlProblem_.finalCost.add(finalCost);

  std::array<double, 6> timeTraj;
  std::array<Eigen::Vector2d, 6> stateTraj, inputTraj;
  for (size_t i = 0; i < 6; ++i) {
    timeTraj[i] = static_cast<double>(i);
    stateTraj[i].setZero();
    inputTraj[i].setZero();
  }
  solver.setDesireTrajectory(timeTraj, stateTraj, inputTraj);

  Eigen::Vector2d initState;
  initState << 1.0, 0.5;

  EXPECT_NO_THROW(solver.run(0.0, initState));
}

// 验证 rollout 统计出的中间代价和终端代价与二次跟踪代价一致。
TEST(iLQREndToEndTest, RolloutMetricsMatchQuadraticTrackingCosts) {
  using Descriptor =
      iLQRDescriptor<double, TranscriptionConfig<Dimensions<2, 2>, Horizon<2>>>;
  using Solver = iLQR<Descriptor>;

  Eigen::Matrix2d A = Eigen::Matrix2d::Zero();
  Eigen::Matrix2d B = Eigen::Matrix2d::Identity();
  LinearSystemDynamics<double, 2, 2> dynamics(A, B);
  DefaultInitializer<double, 2, 2> initializer;
  DDPSettings<double> ddp_setting;
  Solver solver(ddp_setting, &dynamics, &initializer);

  QuadraticStateInputCost<double, 2, 2, 3> runningCost(
      Eigen::Matrix2d::Identity(), Eigen::Matrix2d::Identity());
  QuadraticStateCost<double, 2, 3> finalCost(2.0 * Eigen::Matrix2d::Identity());
  solver.optimalControlProblem_.cost.add(runningCost);
  solver.optimalControlProblem_.finalCost.add(finalCost);

  std::array<double, 3> timeTraj = {0.0, 1.0, 2.0};
  std::array<Eigen::Vector2d, 3> stateTraj, inputTraj;
  for (size_t i = 0; i < 3; ++i) {
    stateTraj[i].setZero();
    inputTraj[i].setZero();
  }
  solver.setDesireTrajectory(timeTraj, stateTraj, inputTraj);

  Solver::PrimalSolution_t primalSolution;
  primalSolution.timeTrajectory_ = timeTraj;
  primalSolution.stateTrajectory_[0] << 1.0, 0.0;
  primalSolution.stateTrajectory_[1] << 2.0, 0.0;
  primalSolution.stateTrajectory_[2] << 3.0, 0.0;
  primalSolution.inputTrajectory_[0] << 1.0, 0.0;
  primalSolution.inputTrajectory_[1] << 0.0, 2.0;
  primalSolution.inputTrajectory_[2].setZero();

  Solver::DualSolution_t dualSolution;
  Solver::ProblemMetrics_t metrics;
  Solver::computeRolloutMetrics(solver.optimalControlProblem_, primalSolution,
                                dualSolution, metrics);

  EXPECT_DOUBLE_EQ(metrics.intermediates[0].cost, 1.0);
  EXPECT_DOUBLE_EQ(metrics.intermediates[1].cost, 4.0);
  EXPECT_DOUBLE_EQ(metrics.final.cost, 9.0);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
