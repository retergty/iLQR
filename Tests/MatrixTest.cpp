#include <gtest/gtest.h>

#include "Matrix.hpp"
#include "Vector.hpp"

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
