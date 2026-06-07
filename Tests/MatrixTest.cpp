#include <gtest/gtest.h>

#include "matrix/Matrix.hpp"
#include "matrix/Vector.hpp"
#include "matrix/Vector3.hpp"

// 验证自实现 Matrix 支持一维 initializer_list 按行优先初始化。
TEST(MatrixTest, SupportsFlatInitializerListConstruction) {
  const matrix::Matrix<double, 2, 3> matrix{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

  EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(matrix(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(matrix(0, 2), 3.0);
  EXPECT_DOUBLE_EQ(matrix(1, 0), 4.0);
  EXPECT_DOUBLE_EQ(matrix(1, 1), 5.0);
  EXPECT_DOUBLE_EQ(matrix(1, 2), 6.0);
}

// 验证自实现 Matrix 支持二维 initializer_list 按行列初始化。
TEST(MatrixTest, SupportsNestedInitializerListConstruction) {
  const matrix::Matrix<double, 2, 2> matrix{{1.0, 2.0}, {3.0, 4.0}};

  EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(matrix(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(matrix(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(matrix(1, 1), 4.0);
}

// 验证已存在对象也可以用一维和二维 initializer_list 重新赋值。
TEST(MatrixTest, SupportsInitializerListAssignment) {
  matrix::Matrix<double, 2, 2> matrix;

  matrix = {5.0, 6.0, 7.0, 8.0};
  EXPECT_DOUBLE_EQ(matrix(0, 0), 5.0);
  EXPECT_DOUBLE_EQ(matrix(0, 1), 6.0);
  EXPECT_DOUBLE_EQ(matrix(1, 0), 7.0);
  EXPECT_DOUBLE_EQ(matrix(1, 1), 8.0);

  matrix = {{9.0, 10.0}, {11.0, 12.0}};
  EXPECT_DOUBLE_EQ(matrix(0, 0), 9.0);
  EXPECT_DOUBLE_EQ(matrix(0, 1), 10.0);
  EXPECT_DOUBLE_EQ(matrix(1, 0), 11.0);
  EXPECT_DOUBLE_EQ(matrix(1, 1), 12.0);
}

// 验证 Matrix 提供 Eigen 风格的静态工厂函数。
TEST(MatrixTest, MatrixSupportsEigenStyleFactories) {
  const matrix::Matrix<double, 2, 2> zero =
      matrix::Matrix<double, 2, 2>::Zero();
  const matrix::Matrix<double, 2, 2> ones =
      matrix::Matrix<double, 2, 2>::Ones();
  const matrix::Matrix<double, 2, 2> constant =
      matrix::Matrix<double, 2, 2>::Constant(3.0);
  const matrix::Matrix<double, 2, 2> identity =
      matrix::Matrix<double, 2, 2>::Identity();

  EXPECT_TRUE(zero.isApprox(
      matrix::Matrix<double, 2, 2>{{0.0, 0.0}, {0.0, 0.0}}, 1e-12));
  EXPECT_TRUE(ones.isApprox(
      matrix::Matrix<double, 2, 2>{{1.0, 1.0}, {1.0, 1.0}}, 1e-12));
  EXPECT_TRUE(constant.isApprox(
      matrix::Matrix<double, 2, 2>{{3.0, 3.0}, {3.0, 3.0}}, 1e-12));
  EXPECT_TRUE(identity.isApprox(
      matrix::Matrix<double, 2, 2>{{1.0, 0.0}, {0.0, 1.0}}, 1e-12));
}

// 验证 Matrix::Random() 生成 Eigen 风格 [-1, 1] 区间随机值。
TEST(MatrixTest, MatrixSupportsRandomFactory) {
  const matrix::Matrix<double, 2, 3> random =
      matrix::Matrix<double, 2, 3>::Random();

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_GE(random(i, j), -1.0);
      EXPECT_LE(random(i, j), 1.0);
    }
  }
}

// 验证 Vector 转发 Matrix 的 initializer_list 构造能力。
TEST(MatrixTest, VectorSupportsInitializerListConstruction) {
  const matrix::Vector<double, 3> vector{1.0, 2.0, 3.0};

  EXPECT_DOUBLE_EQ(vector(0), 1.0);
  EXPECT_DOUBLE_EQ(vector(1), 2.0);
  EXPECT_DOUBLE_EQ(vector(2), 3.0);
}

// 验证 Vector 转发 Matrix 的 initializer_list 赋值能力。
TEST(MatrixTest, VectorSupportsInitializerListAssignment) {
  matrix::Vector<double, 3> vector;

  vector = {4.0, 5.0, 6.0};

  EXPECT_DOUBLE_EQ(vector(0), 4.0);
  EXPECT_DOUBLE_EQ(vector(1), 5.0);
  EXPECT_DOUBLE_EQ(vector(2), 6.0);
}

// 验证固定尺寸 Vector 支持显式跨标量类型构造。
TEST(MatrixTest, VectorSupportsExplicitScalarConversionConstruction) {
  const matrix::Vector<float, 3> source{1.25f, -2.5f, 3.75f};

  const matrix::Vector<double, 3> converted(source);

  EXPECT_TRUE(
      converted.isApprox(matrix::Vector<double, 3>{1.25, -2.5, 3.75}, 1e-12));
}

// 验证 Vector3f 可以显式构造成 Vector3d，兼容 PX4 中的 Vector3d(Vector3f)
// 写法。
TEST(MatrixTest, Vector3SupportsExplicitScalarConversionConstruction) {
  const matrix::Vector3f source{1.25f, -2.5f, 3.75f};

  const matrix::Vector3d converted(source);

  EXPECT_TRUE(converted.isApprox(matrix::Vector3d{1.25, -2.5, 3.75}, 1e-12));
}

// 验证 Vector 提供返回 Vector 类型的静态工厂函数。
TEST(MatrixTest, VectorSupportsEigenStyleFactories) {
  const matrix::Vector<double, 3> zero = matrix::Vector<double, 3>::Zero();
  const matrix::Vector<double, 3> ones = matrix::Vector<double, 3>::Ones();
  const matrix::Vector<double, 3> constant =
      matrix::Vector<double, 3>::Constant(3.0);
  const matrix::Vector<double, 3> identity =
      matrix::Vector<double, 3>::Identity();

  EXPECT_TRUE(zero.isApprox(matrix::Vector<double, 3>{0.0, 0.0, 0.0}, 1e-12));
  EXPECT_TRUE(ones.isApprox(matrix::Vector<double, 3>{1.0, 1.0, 1.0}, 1e-12));
  EXPECT_TRUE(
      constant.isApprox(matrix::Vector<double, 3>{3.0, 3.0, 3.0}, 1e-12));
  EXPECT_TRUE(
      identity.isApprox(matrix::Vector<double, 3>{1.0, 0.0, 0.0}, 1e-12));
}

// 验证 Vector::Random() 返回 Vector 类型并生成 [-1, 1] 区间随机值。
TEST(MatrixTest, VectorSupportsRandomFactory) {
  const matrix::Vector<double, 3> random = matrix::Vector<double, 3>::Random();

  for (int i = 0; i < 3; ++i) {
    EXPECT_GE(random(i), -1.0);
    EXPECT_LE(random(i), 1.0);
  }
}

// 验证 Matrix 支持自定义 Infinity 常量的无穷范数接口。
TEST(MatrixTest, MatrixSupportsInfinityLpNorm) {
  const matrix::Matrix<double, 2, 3> matrix{-1.0, 2.0, -5.0, 4.0, -3.0, 0.5};

  EXPECT_DOUBLE_EQ(matrix.lpNorm<matrix::Infinity>(), 5.0);
}

// 验证 Vector 通过 Matrix 基类支持无穷范数接口。
TEST(MatrixTest, VectorSupportsInfinityLpNorm) {
  const matrix::Vector<double, 3> vector{-1.0, 2.0, -5.0};

  EXPECT_DOUBLE_EQ(vector.lpNorm<matrix::Infinity>(), 5.0);
}

// 验证 Vector 提供 Eigen 风格的固定尺寸 head/tail 包装函数。
TEST(MatrixTest, VectorSupportsHeadAndTail) {
  matrix::Vector<double, 5> vector{1.0, 2.0, 3.0, 4.0, 5.0};

  vector.template head<2>() = matrix::Vector<double, 2>{10.0, 20.0};
  vector.template tail<2>() = matrix::Vector<double, 2>{40.0, 50.0};

  EXPECT_TRUE(vector.isApprox(
      matrix::Vector<double, 5>{10.0, 20.0, 3.0, 40.0, 50.0}, 1e-12));

  const matrix::Vector<double, 5>& constVector = vector;
  const matrix::Vector<double, 2> head = constVector.template head<2>();
  const matrix::Vector<double, 2> tail = constVector.template tail<2>();

  EXPECT_TRUE(head.isApprox(matrix::Vector<double, 2>{10.0, 20.0}, 1e-12));
  EXPECT_TRUE(tail.isApprox(matrix::Vector<double, 2>{40.0, 50.0}, 1e-12));
}

// 验证 Vector 可以直接接收矩阵行/列切片赋值，避免切片到 Vector/Matrix
// 的隐式转换二义性。
TEST(MatrixTest, VectorSupportsRowAndColumnSliceAssignment) {
  matrix::Matrix<double, 3, 2> matrix{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};

  matrix::Vector<double, 3> writableColumn;
  writableColumn = matrix.col(1);
  EXPECT_TRUE(
      writableColumn.isApprox(matrix::Vector<double, 3>{2.0, 4.0, 6.0}, 1e-12));

  const matrix::Matrix<double, 3, 2>& constMatrix = matrix;
  matrix::Vector<double, 3> constColumn;
  constColumn = constMatrix.col(0);
  EXPECT_TRUE(
      constColumn.isApprox(matrix::Vector<double, 3>{1.0, 3.0, 5.0}, 1e-12));

  matrix::Vector<double, 2> writableRow;
  writableRow = matrix.row(2);
  EXPECT_TRUE(writableRow.isApprox(matrix::Vector<double, 2>{5.0, 6.0}, 1e-12));

  matrix::Vector<double, 2> constRow;
  constRow = constMatrix.row(0);
  EXPECT_TRUE(constRow.isApprox(matrix::Vector<double, 2>{1.0, 2.0}, 1e-12));
}

// 验证 Matrix 提供 Eigen 风格的固定尺寸块包装函数。
TEST(MatrixTest, MatrixSupportsFixedSizeBlockWrappers) {
  matrix::Matrix<double, 4, 4> matrix = matrix::Matrix<double, 4, 4>::Zero();

  matrix.template topLeftCorner<2, 2>() =
      matrix::Matrix<double, 2, 2>{{1.0, 2.0}, {3.0, 4.0}};
  matrix.template topRightCorner<2, 2>() =
      matrix::Matrix<double, 2, 2>{{5.0, 6.0}, {7.0, 8.0}};
  matrix.template bottomLeftCorner<2, 2>() =
      matrix::Matrix<double, 2, 2>{{9.0, 10.0}, {11.0, 12.0}};
  matrix.template bottomRightCorner<2, 2>() =
      matrix::Matrix<double, 2, 2>{{13.0, 14.0}, {15.0, 16.0}};
  matrix.template topRows<1>() =
      matrix::Matrix<double, 1, 4>{17.0, 18.0, 19.0, 20.0};
  matrix.template bottomRows<1>() =
      matrix::Matrix<double, 1, 4>{21.0, 22.0, 23.0, 24.0};

  const matrix::Matrix<double, 4, 4> expected{{17.0, 18.0, 19.0, 20.0},
                                              {3.0, 4.0, 7.0, 8.0},
                                              {9.0, 10.0, 13.0, 14.0},
                                              {21.0, 22.0, 23.0, 24.0}};

  EXPECT_TRUE(matrix.isApprox(expected, 1e-12));

  const matrix::Matrix<double, 4, 4>& constMatrix = matrix;
  const matrix::Matrix<double, 2, 2> topLeft =
      constMatrix.template topLeftCorner<2, 2>();
  const matrix::Matrix<double, 2, 2> topRight =
      constMatrix.template topRightCorner<2, 2>();
  const matrix::Matrix<double, 2, 2> bottomLeft =
      constMatrix.template bottomLeftCorner<2, 2>();
  const matrix::Matrix<double, 2, 2> bottomRight =
      constMatrix.template bottomRightCorner<2, 2>();

  EXPECT_TRUE(topLeft.isApprox(
      matrix::Matrix<double, 2, 2>{{17.0, 18.0}, {3.0, 4.0}}, 1e-12));
  EXPECT_TRUE(topRight.isApprox(
      matrix::Matrix<double, 2, 2>{{19.0, 20.0}, {7.0, 8.0}}, 1e-12));
  EXPECT_TRUE(bottomLeft.isApprox(
      matrix::Matrix<double, 2, 2>{{9.0, 10.0}, {21.0, 22.0}}, 1e-12));
  EXPECT_TRUE(bottomRight.isApprox(
      matrix::Matrix<double, 2, 2>{{13.0, 14.0}, {23.0, 24.0}}, 1e-12));
}

// 验证可写切片支持 Eigen 风格 setIdentity()。
TEST(MatrixTest, SliceSupportsSetIdentity) {
  matrix::Matrix<double, 4, 4> matrix = matrix::Matrix<double, 4, 4>::Ones();

  matrix.template topLeftCorner<3, 2>().setIdentity();

  const matrix::Matrix<double, 4, 4> expected{{1.0, 0.0, 1.0, 1.0},
                                              {0.0, 1.0, 1.0, 1.0},
                                              {0.0, 0.0, 1.0, 1.0},
                                              {1.0, 1.0, 1.0, 1.0}};

  EXPECT_TRUE(matrix.isApprox(expected, 1e-12));
}

// 验证切片可直接与矩阵相乘，不需要先构造临时矩阵。
TEST(MatrixTest, SliceSupportsMatrixMultiplication) {
  const matrix::Matrix<double, 4, 3> matrix{
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}, {10.0, 11.0, 12.0}};
  const matrix::Vector<double, 3> vector{1.0, 2.0, 3.0};

  const matrix::Matrix<double, 2, 1> product =
      matrix.template topRows<2>() * vector;
  const matrix::Vector<double, 2> expected{14.0, 32.0};

  EXPECT_TRUE(product.isApprox(expected, 1e-12));
}

// 验证切片的非原地加减和标量乘法不依赖临时切片矩阵。
TEST(MatrixTest, SliceSupportsNonMutatingArithmetic) {
  matrix::Matrix<double, 2, 4> matrix{{1.0, 2.0, 5.0, 6.0},
                                      {3.0, 4.0, 7.0, 8.0}};
  const matrix::Matrix<double, 2, 4>& constMatrix = matrix;
  const matrix::Matrix<double, 2, 2> rhs{{10.0, 20.0}, {30.0, 40.0}};

  const matrix::Matrix<double, 2, 2> plusMatrix =
      matrix.template topLeftCorner<2, 2>() + rhs;
  const matrix::Matrix<double, 2, 2> minusMatrix =
      matrix.template topLeftCorner<2, 2>() - rhs;
  const matrix::Matrix<double, 2, 2> plusSlice =
      matrix.template topLeftCorner<2, 2>() +
      constMatrix.template topRightCorner<2, 2>();
  const matrix::Matrix<double, 2, 2> minusSlice =
      matrix.template topRightCorner<2, 2>() -
      constMatrix.template topLeftCorner<2, 2>();
  const matrix::Matrix<double, 2, 2> plusScalar =
      matrix.template topLeftCorner<2, 2>() + 1.0;
  const matrix::Matrix<double, 2, 2> minusScalar =
      matrix.template topLeftCorner<2, 2>() - 1.0;
  const matrix::Matrix<double, 2, 2> scaled =
      matrix.template topLeftCorner<2, 2>() * 2.0;

  EXPECT_TRUE(plusMatrix.isApprox(
      matrix::Matrix<double, 2, 2>{{11.0, 22.0}, {33.0, 44.0}}, 1e-12));
  EXPECT_TRUE(minusMatrix.isApprox(
      matrix::Matrix<double, 2, 2>{{-9.0, -18.0}, {-27.0, -36.0}}, 1e-12));
  EXPECT_TRUE(plusSlice.isApprox(
      matrix::Matrix<double, 2, 2>{{6.0, 8.0}, {10.0, 12.0}}, 1e-12));
  EXPECT_TRUE(minusSlice.isApprox(
      matrix::Matrix<double, 2, 2>{{4.0, 4.0}, {4.0, 4.0}}, 1e-12));
  EXPECT_TRUE(plusScalar.isApprox(
      matrix::Matrix<double, 2, 2>{{2.0, 3.0}, {4.0, 5.0}}, 1e-12));
  EXPECT_TRUE(minusScalar.isApprox(
      matrix::Matrix<double, 2, 2>{{0.0, 1.0}, {2.0, 3.0}}, 1e-12));
  EXPECT_TRUE(scaled.isApprox(
      matrix::Matrix<double, 2, 2>{{2.0, 4.0}, {6.0, 8.0}}, 1e-12));
}

// 验证 Vector 不隐藏 Matrix 的矩阵乘法重载。
TEST(MatrixTest, VectorSupportsMatrixMultiplication) {
  const matrix::Vector<double, 2> vector{2.0, 3.0};
  const matrix::Matrix<double, 1, 3> row{4.0, 5.0, 6.0};

  const matrix::Matrix<double, 2, 3> product = vector * row;
  const matrix::Matrix<double, 2, 3> expected{{8.0, 10.0, 12.0},
                                              {12.0, 15.0, 18.0}};

  EXPECT_TRUE(product.isApprox(expected, 1e-12));
}

// 验证 Matrix 可以直接计算 A^T * B，避免显式构造转置矩阵。
TEST(MatrixTest, MatrixSupportsTransposeMultiply) {
  const matrix::Matrix<double, 3, 2> lhs{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  const matrix::Matrix<double, 3, 2> rhs{
      {7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}};

  const matrix::Matrix<double, 2, 2> product = lhs.transposeMultiply(rhs);
  const matrix::Matrix<double, 2, 2> expected = lhs.transpose() * rhs;

  EXPECT_TRUE(product.isApprox(expected, 1e-12));
}

// 验证 Matrix 可以直接累加 A^T * B，不产生乘积临时矩阵。
TEST(MatrixTest, MatrixSupportsAddTransposeProduct) {
  const matrix::Matrix<double, 3, 2> lhs{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  const matrix::Matrix<double, 3, 1> rhs{7.0, 8.0, 9.0};
  matrix::Matrix<double, 2, 1> accumulated{10.0, 20.0};

  const matrix::Matrix<double, 2, 1> expected =
      accumulated + lhs.transpose() * rhs;
  accumulated.addTransposeProduct(lhs, rhs);

  EXPECT_TRUE(accumulated.isApprox(expected, 1e-12));
}

// 验证 Vector * Vector 保持点积语义。
TEST(MatrixTest, VectorMultiplicationWithVectorKeepsDotProduct) {
  const matrix::Vector<double, 2> lhs{2.0, 3.0};
  const matrix::Vector<double, 2> rhs{4.0, 5.0};

  EXPECT_DOUBLE_EQ(lhs * rhs, 23.0);
}

// 验证 1 维 Vector 与 1x1 Matrix 相乘仍走矩阵乘法。
TEST(MatrixTest, OneDimensionalVectorSupportsMatrixMultiplication) {
  const matrix::Vector<double, 1> vector{2.0};
  const matrix::Matrix<double, 1, 1> matrix{3.0};

  const matrix::Matrix<double, 1, 1> product = vector * matrix;

  EXPECT_TRUE(product.isApprox(matrix::Matrix<double, 1, 1>{6.0}, 1e-12));
}

// 验证 Matrix、Vector 和 Slice 支持 Eigen 风格 squaredNorm()。
TEST(MatrixTest, SupportsSquaredNorm) {
  matrix::Matrix<double, 2, 2> matrix{{1.0, -2.0}, {3.0, -4.0}};
  const matrix::Vector<double, 3> vector{1.0, -2.0, 3.0};

  EXPECT_DOUBLE_EQ(matrix.squaredNorm(), 30.0);
  EXPECT_DOUBLE_EQ(vector.squaredNorm(), 14.0);
  EXPECT_DOUBLE_EQ((matrix.template topLeftCorner<1, 2>().squaredNorm()), 5.0);
}

// 验证 Matrix 支持 Eigen 风格 Frobenius norm()。
TEST(MatrixTest, MatrixSupportsNorm) {
  const matrix::Matrix<double, 2, 2> matrix{{1.0, -2.0}, {3.0, -4.0}};
  const matrix::Vector<double, 2> lhs{3.0, 4.0};
  const matrix::Vector<double, 2> rhs{1.0, 1.0};

  EXPECT_DOUBLE_EQ(matrix.norm(), std::sqrt(30.0));
  EXPECT_DOUBLE_EQ((lhs - rhs).norm(), std::sqrt(13.0));
}

// 验证 Matrix、Vector 和 Slice 支持 Eigen 风格 isZero()。
TEST(MatrixTest, SupportsIsZero) {
  matrix::Matrix<double, 2, 2> matrix = matrix::Matrix<double, 2, 2>::Zero();
  matrix(0, 1) = 1e-13;

  EXPECT_TRUE(matrix.isZero(1e-12));
  EXPECT_FALSE(matrix.isZero(1e-14));

  const matrix::Vector<double, 3> vector{0.0, -1e-13, 0.0};
  EXPECT_TRUE(vector.isZero(1e-12));
  EXPECT_FALSE(vector.isZero(1e-14));

  EXPECT_TRUE((matrix.template topLeftCorner<1, 1>().isZero(1e-12)));
}

// 验证 Matrix 和 Vector 支持 Eigen 风格对象级 swap()。
TEST(MatrixTest, SupportsObjectSwap) {
  matrix::Matrix<double, 2, 2> lhs{{1.0, 2.0}, {3.0, 4.0}};
  matrix::Matrix<double, 2, 2> rhs{{5.0, 6.0}, {7.0, 8.0}};

  lhs.swap(rhs);

  EXPECT_TRUE(lhs.isApprox(matrix::Matrix<double, 2, 2>{{5.0, 6.0}, {7.0, 8.0}},
                           1e-12));
  EXPECT_TRUE(rhs.isApprox(matrix::Matrix<double, 2, 2>{{1.0, 2.0}, {3.0, 4.0}},
                           1e-12));

  matrix::Vector<double, 3> lhsVector{1.0, 2.0, 3.0};
  matrix::Vector<double, 3> rhsVector{4.0, 5.0, 6.0};

  lhsVector.swap(rhsVector);

  EXPECT_TRUE(
      lhsVector.isApprox(matrix::Vector<double, 3>{4.0, 5.0, 6.0}, 1e-12));
  EXPECT_TRUE(
      rhsVector.isApprox(matrix::Vector<double, 3>{1.0, 2.0, 3.0}, 1e-12));
}

// 验证自实现 Matrix 提供 Eigen 风格的 isApprox 成员接口。
TEST(MatrixTest, SupportsIsApprox) {
  const matrix::Matrix<double, 2, 2> lhs{{1.0, 2.0}, {3.0, 4.0}};
  const matrix::Matrix<double, 2, 2> close{{1.0 + 1e-13, 2.0},
                                           {3.0, 4.0 - 1e-13}};
  const matrix::Matrix<double, 2, 2> far{{1.0 + 1e-6, 2.0}, {3.0, 4.0}};

  EXPECT_TRUE(lhs.isApprox(close, 1e-12));
  EXPECT_FALSE(lhs.isApprox(far, 1e-12));
}

// 验证 Vector 继承 Matrix 的 isApprox 接口，便于替换 Eigen::Vector 用法。
TEST(MatrixTest, VectorSupportsIsApprox) {
  const matrix::Vector<double, 3> lhs(
      matrix::Matrix<double, 3, 1>{1.0, 2.0, 3.0});
  const matrix::Vector<double, 3> rhs(
      matrix::Matrix<double, 3, 1>{1.0, 2.0 + 1e-13, 3.0});

  EXPECT_TRUE(lhs.isApprox(rhs, 1e-12));
}
