#include <gtest/gtest.h>

#include <Eigen/Cholesky>
#include <array>

#include "Controller/LinearController.hpp"
#include "Initialization/DefaultInitializer.hpp"
#include "ModelData/Multiplier.hpp"
#include "Tests/Include/MatrixEigenConversion.hpp"
#include "iLQR/HessianCorrection.hpp"
#include "iLQR/iLQRDescriptor.hpp"
#include "matrix/CholeskyDecomposition.hpp"
#include "matrix/Types.hpp"

using test_tools::matrix_eigen_conversion::toEigenMatrix;

// 验证对角 Hessian 修正只会平移对角项。
TEST(UtilityTest, ShiftHessianAddsOnlyDiagonalShift) {
  Matrix<double, 2, 2> matrix{{1.0, 2.0}, {3.0, 4.0}};

  shiftHessian(HessianCorrectionStrategy::DIAGONAL_SHIFT, matrix, 0.25);

  const Matrix<double, 2, 2> expected{{1.25, 2.0}, {3.0, 4.25}};
  EXPECT_TRUE(matrix.isApprox(expected, 1e-12));
}

// 验证自实现 Cholesky 分解支持矩阵右端项求解。
TEST(UtilityTest, CholeskyDecompositionSolvesMatrixRightHandSide) {
  Matrix<double, 3, 3> matrix{
      {6.0, 2.0, 1.0}, {2.0, 5.0, 2.0}, {1.0, 2.0, 4.0}};
  Matrix<double, 3, 2> rhs{{1.0, 2.0}, {-3.0, 4.0}, {5.0, -6.0}};

  CholeskyDecomposition<double, 3> cholesky;
  ASSERT_TRUE(cholesky.Decomposition(matrix));

  Matrix<double, 3, 2> solution;
  cholesky.Solve(solution, rhs);

  const Eigen::Matrix<double, 3, 2> expected =
      toEigenMatrix(matrix).llt().solve(toEigenMatrix(rhs));
  EXPECT_TRUE(toEigenMatrix(solution).isApprox(expected, 1e-12));
}

// 验证 Cholesky 只存下三角时仍能重构 L、L^T 并通过 const 对象求解。
TEST(UtilityTest, CholeskyDecompositionAccessorsAreConstCorrect) {
  Matrix<double, 3, 3> matrix{
      {6.0, 2.0, 1.0}, {2.0, 5.0, 2.0}, {1.0, 2.0, 4.0}};
  Vector<double, 3> rhs{1.0, -3.0, 5.0};

  CholeskyDecomposition<double, 3> cholesky;
  ASSERT_TRUE(cholesky.Decomposition(matrix));
  const CholeskyDecomposition<double, 3>& constCholesky = cholesky;

  const Matrix<double, 3, 3> l = constCholesky.GetMatrixL();
  const Matrix<double, 3, 3> lt = constCholesky.GetMatrixLT();
  EXPECT_TRUE((l * lt).isApprox(matrix, 1e-12));

  Vector<double, 3> solution;
  constCholesky.Solve(solution, rhs);
  const Eigen::Matrix<double, 3, 1> expected =
      toEigenMatrix(matrix).llt().solve(toEigenMatrix(rhs));
  EXPECT_TRUE(toEigenMatrix(solution).isApprox(expected, 1e-12));
}

// 验证 LinearController 在计算输入前会同时插值 bias 和 gain。
TEST(UtilityTest, LinearControllerInterpolatesBiasAndGain) {
  LinearController<double, 2, 1, 2> controller;
  controller.timeStamp_[0] = 0.0;
  controller.timeStamp_[1] = 2.0;
  controller.biasArray_[0] = {1.0};
  controller.biasArray_[1] = {3.0};
  controller.gainArray_[0] = {1.0, 0.0};
  controller.gainArray_[1] = {3.0, 2.0};

  const Vector<double, 2> state{2.0, 4.0};

  const Vector<double, 1> expected{10.0};

  EXPECT_TRUE(controller.computeInput(1.0, state).isApprox(expected, 1e-12));
}

// 验证 clear() 和 swap() 后控制器的空与非空状态符合预期。
TEST(UtilityTest, LinearControllerClearAndSwap) {
  LinearController<double, 2, 1, 2> first;
  first.timeStamp_[0] = 0.0;
  first.timeStamp_[1] = 1.0;
  first.biasArray_[0] = {10.0};
  first.biasArray_[1] = {20.0};
  first.gainArray_[0] = Matrix<double, 1, 2>::Zero();
  first.gainArray_[1] = Matrix<double, 1, 2>::Zero();

  LinearController<double, 2, 1, 2> second;
  second.clear();

  EXPECT_FALSE(first.empty());
  EXPECT_TRUE(second.empty());

  first.swap(second);
  EXPECT_TRUE(first.empty());
  EXPECT_FALSE(second.empty());

  second.clear();
  EXPECT_TRUE(second.empty());
}

// 验证单个乘子插值会同时插值 penalty 和 lagrangian 字段。
TEST(UtilityTest, MultiplierInterpolationInterpolatesAllFields) {
  using Multiplier_t = Multiplier<double, 1>;
  std::array<Multiplier_t, 3> multipliers = {
      Multiplier_t{1.0, Vector<double, 1>{10.0}},
      Multiplier_t{3.0, Vector<double, 1>{30.0}},
      Multiplier_t{5.0, Vector<double, 1>{50.0}},
  };

  const auto interpolated =
      LinearInterpolation::interpolate<double, 1, 3>({0, 0.25}, multipliers);

  EXPECT_DOUBLE_EQ(interpolated.penalty, 2.5);
  EXPECT_DOUBLE_EQ(interpolated.lagrangian(0), 25.0);
}

// 验证乘子集合插值会分别处理各类约束乘子。
TEST(UtilityTest, MultiplierCollectionInterpolationInterpolatesEachCategory) {
  using OneDimGroup = ConstraintGroupLayout<ConstraintTerm<1>>;
  using Layout =
      StageConstraintLayout<OneDimGroup, OneDimGroup, OneDimGroup, OneDimGroup>;
  using Collection = MultiplierCollection<double, Layout>;
  std::array<Collection, 2> trajectory;

  std::get<0>(trajectory[0].stateEq.terms) = {1.0, Vector<double, 1>{10.0}};
  std::get<0>(trajectory[0].stateIneq.terms) = {2.0, Vector<double, 1>{20.0}};
  std::get<0>(trajectory[0].stateInputEq.terms) = {3.0,
                                                   Vector<double, 1>{30.0}};
  std::get<0>(trajectory[0].stateInputIneq.terms) = {4.0,
                                                     Vector<double, 1>{40.0}};
  std::get<0>(trajectory[1].stateEq.terms) = {11.0, Vector<double, 1>{110.0}};
  std::get<0>(trajectory[1].stateIneq.terms) = {12.0, Vector<double, 1>{120.0}};
  std::get<0>(trajectory[1].stateInputEq.terms) = {13.0,
                                                   Vector<double, 1>{130.0}};
  std::get<0>(trajectory[1].stateInputIneq.terms) = {14.0,
                                                     Vector<double, 1>{140.0}};

  const Collection interpolated =
      LinearInterpolation::interpolate<double, Layout, 2>({0, 0.25},
                                                          trajectory);

  EXPECT_DOUBLE_EQ(std::get<0>(interpolated.stateEq.terms).penalty, 8.5);
  EXPECT_DOUBLE_EQ(std::get<0>(interpolated.stateEq.terms).lagrangian(0), 85.0);
  EXPECT_DOUBLE_EQ(std::get<0>(interpolated.stateIneq.terms).penalty, 9.5);
  EXPECT_DOUBLE_EQ(std::get<0>(interpolated.stateIneq.terms).lagrangian(0),
                   95.0);
  EXPECT_DOUBLE_EQ(std::get<0>(interpolated.stateInputEq.terms).penalty, 10.5);
  EXPECT_DOUBLE_EQ(std::get<0>(interpolated.stateInputEq.terms).lagrangian(0),
                   105.0);
  EXPECT_DOUBLE_EQ(std::get<0>(interpolated.stateInputIneq.terms).penalty,
                   11.5);
  EXPECT_DOUBLE_EQ(std::get<0>(interpolated.stateInputIneq.terms).lagrangian(0),
                   115.0);
}

// 验证默认初始化器保持状态不变并把输入初值清零。
TEST(UtilityTest, DefaultInitializerKeepsStateAndZerosInput) {
  DefaultInitializer<double, 2, 1> initializer;
  const Vector<double, 2> state{-1.0, 2.0};
  Vector<double, 1> input{9.0};
  Vector<double, 2> nextState = Vector<double, 2>::Zero();

  initializer.compute(0.0, state, 1.0, input, nextState);

  EXPECT_DOUBLE_EQ(input(0), 0.0);
  EXPECT_TRUE(nextState.isApprox(state, 1e-12));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
