#include <gtest/gtest.h>

#include <array>

#include "ChangeOfInputVariables.hpp"
#include "LinearApproximation.hpp"
#include "LinearQuadraticApproximator.hpp"
#include "LinearSystemDynamics.hpp"
#include "QuadraticApproximation.hpp"
#include "QuadraticStateCost.hpp"
#include "iLQRDescriptor.hpp"

// 验证输入变量变换只影响二次近似里与输入相关的项。
TEST(ApproximationTest,
     ChangeOfInputVariablesTransformsQuadraticApproximation) {
  ScalarFunctionQuadraticApproximation<double, 2, 2> approximation;
  approximation.dfdxx << 13.0, 17.0, 19.0, 23.0;
  approximation.dfdux << 1.0, 2.0, 3.0, 4.0;
  approximation.dfduu << 5.0, 1.0, 1.0, 3.0;
  approximation.dfdx << 29.0, 31.0;
  approximation.dfdu << 7.0, 11.0;
  approximation.f = 37.0;

  const auto originalDfdxx = approximation.dfdxx;
  const auto originalDfdux = approximation.dfdux;
  const auto originalDfduu = approximation.dfduu;
  const auto originalDfdx = approximation.dfdx;
  const auto originalDfdu = approximation.dfdu;
  const double originalF = approximation.f;

  Eigen::Matrix2d Pu;
  Pu << 2.0, -1.0, 0.0, 3.0;

  changeOfInputVariables(approximation, Pu);

  EXPECT_TRUE(
      approximation.dfdux.isApprox(Pu.transpose() * originalDfdux, 1e-12));
  EXPECT_TRUE(
      approximation.dfduu.isApprox(Pu.transpose() * originalDfduu * Pu, 1e-12));
  EXPECT_TRUE(
      approximation.dfdu.isApprox(Pu.transpose() * originalDfdu, 1e-12));
  EXPECT_TRUE(approximation.dfdxx.isApprox(originalDfdxx, 1e-12));
  EXPECT_TRUE(approximation.dfdx.isApprox(originalDfdx, 1e-12));
  EXPECT_DOUBLE_EQ(approximation.f, originalF);
}

// 验证输入变量变换只影响线性近似里的输入 Jacobian。
TEST(ApproximationTest, ChangeOfInputVariablesTransformsLinearApproximation) {
  VectorFunctionLinearApproximation<double, 2, 2, 2> approximation;
  approximation.dfdx << 1.0, 2.0, 3.0, 4.0;
  approximation.dfdu << 5.0, 6.0, 7.0, 8.0;
  approximation.f << 9.0, 10.0;

  const auto originalDfdx = approximation.dfdx;
  const auto originalDfdu = approximation.dfdu;
  const auto originalF = approximation.f;

  Eigen::Matrix2d Pu;
  Pu << 2.0, -1.0, 0.0, 3.0;

  changeOfInputVariables(approximation, Pu);

  EXPECT_TRUE(approximation.dfdu.isApprox(originalDfdu * Pu, 1e-12));
  EXPECT_TRUE(approximation.dfdx.isApprox(originalDfdx, 1e-12));
  EXPECT_TRUE(approximation.f.isApprox(originalF, 1e-12));
}

// 验证线性近似和二次近似结构体的基础代数运算。
TEST(ApproximationTest, LinearAndQuadraticApproximationUtilities) {
  auto linear = ScalarFunctionLinearApproximation<double, 2, 1>::Zero();
  EXPECT_TRUE(linear.dfdx.isZero(1e-12));
  EXPECT_TRUE(linear.dfdu.isZero(1e-12));
  EXPECT_DOUBLE_EQ(linear.f, 0.0);

  ScalarFunctionLinearApproximation<double, 2, 1> linearIncrement;
  linearIncrement.dfdx << 1.0, 2.0;
  linearIncrement.dfdu << 3.0;
  linearIncrement.f = 4.0;
  linear += linearIncrement;
  linear *= 2.0;

  Eigen::Vector2d expectedLinearDfdx;
  expectedLinearDfdx << 2.0, 4.0;
  Eigen::Matrix<double, 1, 1> expectedLinearDfdu;
  expectedLinearDfdu << 6.0;
  EXPECT_TRUE(linear.dfdx.isApprox(expectedLinearDfdx, 1e-12));
  EXPECT_TRUE(linear.dfdu.isApprox(expectedLinearDfdu, 1e-12));
  EXPECT_DOUBLE_EQ(linear.f, 8.0);

  ScalarFunctionQuadraticApproximation<double, 2, 1> quadratic;
  quadratic.setZero();
  quadratic.dfdxx << 1.0, 2.0, 3.0, 4.0;
  quadratic.dfdux << 5.0, 6.0;
  quadratic.dfduu << 7.0;
  quadratic.dfdx << 8.0, 9.0;
  quadratic.dfdu << 10.0;
  quadratic.f = 11.0;

  const auto scaled = 0.5 * quadratic;

  EXPECT_TRUE(scaled.dfdxx.isApprox(0.5 * quadratic.dfdxx, 1e-12));
  EXPECT_TRUE(scaled.dfdux.isApprox(0.5 * quadratic.dfdux, 1e-12));
  EXPECT_TRUE(scaled.dfduu.isApprox(0.5 * quadratic.dfduu, 1e-12));
  EXPECT_TRUE(scaled.dfdx.isApprox(0.5 * quadratic.dfdx, 1e-12));
  EXPECT_TRUE(scaled.dfdu.isApprox(0.5 * quadratic.dfdu, 1e-12));
  EXPECT_DOUBLE_EQ(scaled.f, 5.5);
}

// 验证 LQ 近似器对中间时刻代价和动力学线性化的计算结果。
TEST(ApproximationTest,
     LinearQuadraticApproximatorComputesIntermediateCostAndLQ) {
  using Descriptor =
      iLQRDescriptor<double, TranscriptionConfig<Dimensions<2, 2>, Horizon<2>>>;
  using Approximator = LinearQuadraticApproximator<Descriptor>;

  Eigen::Matrix2d A;
  A << 1.0, 2.0, 3.0, 4.0;
  Eigen::Matrix2d B = Eigen::Matrix2d::Identity();
  LinearSystemDynamics<double, 2, 2> dynamics(A, B);

  Approximator::OptimalControlProblem_t problem;
  problem.dynamicsPtr = &dynamics;

  Eigen::Matrix2d Q;
  Q << 2.0, 0.0, 0.0, 4.0;
  Eigen::Matrix2d R;
  R << 3.0, 0.0, 0.0, 5.0;
  Eigen::Matrix2d QState;
  QState << 7.0, 0.0, 0.0, 11.0;
  QuadraticStateInputCost<double, 2, 2, 3> stateInputCost(Q, R);
  QuadraticStateCost<double, 2, 3> stateCost(QState);
  problem.cost.add(stateInputCost);
  problem.stateCost.add(stateCost);

  Approximator::TargetTrajectories_t targetTrajectory;
  targetTrajectory.timeTrajectory = {0.0, 1.0, 2.0};
  targetTrajectory.stateTrajectory[0] << 0.0, 0.0;
  targetTrajectory.stateTrajectory[1] << 1.0, 1.0;
  targetTrajectory.stateTrajectory[2] << 2.0, 2.0;
  targetTrajectory.inputTrajectory[0] << 0.0, 0.0;
  targetTrajectory.inputTrajectory[1] << 1.0, -1.0;
  targetTrajectory.inputTrajectory[2] << 2.0, -2.0;

  Eigen::Vector2d state;
  state << 2.0, 3.0;
  Eigen::Vector2d input;
  input << 4.0, -1.0;

  const double cost =
      Approximator::computeCost(problem, targetTrajectory, 0.5, state, input);
  const auto costApproximation = Approximator::approximateCost(
      problem, targetTrajectory, 0.5, state, input);
  Approximator::IntermediateMultiplierCollection_t multipliers;
  const auto modelData = Approximator::approximateIntermediateLQ(
      problem, targetTrajectory, 0.5, state, input, multipliers);

  Eigen::Vector2d expectedDfdx;
  expectedDfdx << 13.5, 37.5;
  Eigen::Vector2d expectedDfdu;
  expectedDfdu << 10.5, -2.5;
  Eigen::Matrix2d expectedDfdxx;
  expectedDfdxx << 9.0, 0.0, 0.0, 15.0;
  Eigen::Vector2d expectedDynamics;
  expectedDynamics << 12.0, 17.0;

  EXPECT_DOUBLE_EQ(cost, 76.0);
  EXPECT_DOUBLE_EQ(costApproximation.f, 76.0);
  EXPECT_TRUE(costApproximation.dfdx.isApprox(expectedDfdx, 1e-12));
  EXPECT_TRUE(costApproximation.dfdu.isApprox(expectedDfdu, 1e-12));
  EXPECT_TRUE(costApproximation.dfdxx.isApprox(expectedDfdxx, 1e-12));
  EXPECT_TRUE(costApproximation.dfduu.isApprox(R, 1e-12));
  EXPECT_TRUE(costApproximation.dfdux.isZero(1e-12));

  EXPECT_TRUE(modelData.dynamics.f.isApprox(expectedDynamics, 1e-12));
  EXPECT_TRUE(modelData.dynamics.dfdx.isApprox(A, 1e-12));
  EXPECT_TRUE(modelData.dynamics.dfdu.isApprox(B, 1e-12));
  EXPECT_TRUE(modelData.cost.dfdx.isApprox(expectedDfdx, 1e-12));
  EXPECT_TRUE(modelData.cost.dfdu.isApprox(expectedDfdu, 1e-12));
  EXPECT_TRUE(modelData.cost.dfdxx.isApprox(expectedDfdxx, 1e-12));
  EXPECT_TRUE(modelData.cost.dfduu.isApprox(R, 1e-12));
}

// 验证 LQ 近似器对终端代价近似和终端模型数据的计算结果。
TEST(ApproximationTest, LinearQuadraticApproximatorComputesFinalCostAndLQ) {
  using Descriptor =
      iLQRDescriptor<double, TranscriptionConfig<Dimensions<2, 2>, Horizon<2>>>;
  using Approximator = LinearQuadraticApproximator<Descriptor>;

  Approximator::OptimalControlProblem_t problem;
  Eigen::Matrix2d QFinal;
  QFinal << 2.0, 0.0, 0.0, 6.0;
  QuadraticStateCost<double, 2, 3> finalCost(QFinal);
  problem.finalCost.add(finalCost);

  Approximator::TargetTrajectories_t targetTrajectory;
  targetTrajectory.timeTrajectory = {0.0, 1.0, 2.0};
  targetTrajectory.stateTrajectory[0] << 0.0, 0.0;
  targetTrajectory.stateTrajectory[1] << 1.0, 1.0;
  targetTrajectory.stateTrajectory[2] << 2.0, 2.0;

  Eigen::Vector2d state;
  state << 3.5, -0.5;

  Approximator::FinalMultiplierCollection_t multipliers;
  const double cost =
      Approximator::computeFinalCost(problem, targetTrajectory, 1.5, state);
  const auto costApproximation =
      Approximator::approximateFinalCost(problem, targetTrajectory, 1.5, state);
  const auto modelData = Approximator::approximateFinalLQ(
      problem, targetTrajectory, 1.5, state, multipliers);

  Eigen::Vector2d expectedDfdx;
  expectedDfdx << 4.0, -12.0;

  EXPECT_DOUBLE_EQ(cost, 16.0);
  EXPECT_DOUBLE_EQ(costApproximation.f, 16.0);
  EXPECT_TRUE(costApproximation.dfdx.isApprox(expectedDfdx, 1e-12));
  EXPECT_TRUE(costApproximation.dfdxx.isApprox(QFinal, 1e-12));

  EXPECT_TRUE(modelData.dynamics.f.isZero(1e-12));
  EXPECT_TRUE(modelData.dynamics.dfdx.isZero(1e-12));
  EXPECT_TRUE(modelData.dynamics.dfdu.isZero(1e-12));
  EXPECT_DOUBLE_EQ(modelData.cost.f, 16.0);
  EXPECT_TRUE(modelData.cost.dfdx.isApprox(expectedDfdx, 1e-12));
  EXPECT_TRUE(modelData.cost.dfdxx.isApprox(QFinal, 1e-12));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
