/**
 * @file LinearInterpolationTest.cpp
 * @brief 线性插值工具测试：固定步长 timeSegment 与 interpolate。
 */
#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "Matrix/Types.hpp"
#include "Misc/LinearInterpolation.hpp"

using namespace LinearInterpolation;

// 验证固定步长 timeSegment 返回的左端索引和插值权重正确。
TEST(LinearInterpolationTest, TimeSegment) {
  std::array<double, 4> times = {0.0, 1.0, 2.0, 3.0};

  auto [idx0, alpha0] = timeSegment(0.5, times);
  EXPECT_EQ(idx0, 0);
  EXPECT_DOUBLE_EQ(alpha0,
                   0.5);  // alpha = (t1 - t) / (t1 - t0) = (1 - 0.5)/1 = 0.5。

  auto [idx1, alpha1] = timeSegment(1.5, times);
  EXPECT_EQ(idx1, 1);
  EXPECT_DOUBLE_EQ(alpha1, 0.5);

  auto [idx2, alpha2] = timeSegment(0.0, times);
  EXPECT_EQ(idx2, 0);
  EXPECT_DOUBLE_EQ(alpha2, 1.0);

  auto [idx3, alpha3] = timeSegment(3.0, times);
  EXPECT_EQ(idx3, 2);
  EXPECT_DOUBLE_EQ(alpha3, 0.0);
}

// 验证查询时间落在采样节点时，固定步长算法选择该节点作为新区间左端。
TEST(LinearInterpolationTest, TimeSegmentAtInteriorKnot) {
  std::array<double, 4> times = {0.0, 1.0, 2.0, 3.0};

  auto [idx, alpha] = timeSegment(1.0, times);
  EXPECT_EQ(idx, 1);
  EXPECT_DOUBLE_EQ(alpha, 1.0);
}

// 验证超出时间范围时固定步长 timeSegment 会钳制到边界区间。
TEST(LinearInterpolationTest, TimeSegmentOutOfRangeClampsToBoundary) {
  std::array<double, 4> times = {0.0, 1.0, 2.0, 3.0};

  auto [idx0, alpha0] = timeSegment(-10.0, times);
  EXPECT_EQ(idx0, 0);
  EXPECT_DOUBLE_EQ(alpha0, 1.0);

  auto [idx1, alpha1] = timeSegment(10.0, times);
  EXPECT_EQ(idx1, 2);
  EXPECT_DOUBLE_EQ(alpha1, 0.0);
}

// 验证向量样本在线性插值节点和中点处的结果正确。
TEST(LinearInterpolationTest, InterpolateVector) {
  std::array<double, 3> times = {0.0, 1.0, 2.0};
  std::array<Vector<double, 2>, 3> data = {Vector<double, 2>{0.0, 0.0},
                                           Vector<double, 2>{1.0, 2.0},
                                           Vector<double, 2>{2.0, 4.0}};

  Vector<double, 2> r0 = interpolate(0.0, times, data);
  EXPECT_TRUE(r0.isApprox(data[0], 1e-10));

  Vector<double, 2> r1 = interpolate(1.0, times, data);
  EXPECT_TRUE(r1.isApprox(data[1], 1e-10));

  Vector<double, 2> r05 = interpolate(0.5, times, data);
  const Vector<double, 2> expected05{0.5, 1.0};
  EXPECT_TRUE(r05.isApprox(expected05, 1e-10));

  Vector<double, 2> r15 = interpolate(1.5, times, data);
  const Vector<double, 2> expected15{1.5, 3.0};
  EXPECT_TRUE(r15.isApprox(expected15, 1e-10));
}

// 验证超出采样范围时插值结果会钳制到边界值。
TEST(LinearInterpolationTest, InterpolateOutOfRangeClampsToBoundary) {
  std::array<double, 3> times = {0.0, 1.0, 2.0};
  std::array<Vector<double, 2>, 3> data = {Vector<double, 2>{-1.0, 2.0},
                                           Vector<double, 2>{0.0, 4.0},
                                           Vector<double, 2>{3.0, 8.0}};

  EXPECT_TRUE(interpolate(-10.0, times, data).isApprox(data[0], 1e-10));
  EXPECT_TRUE(interpolate(10.0, times, data).isApprox(data[2], 1e-10));
}

// 验证线性插值同样适用于矩阵样本。
TEST(LinearInterpolationTest, InterpolateMatrix) {
  std::array<double, 2> times = {0.0, 2.0};
  std::array<Matrix<double, 2, 2>, 2> data = {
      Matrix<double, 2, 2>{{1.0, 2.0}, {3.0, 4.0}},
      Matrix<double, 2, 2>{{5.0, 6.0}, {7.0, 8.0}}};

  const Matrix<double, 2, 2> expected{{3.0, 4.0}, {5.0, 6.0}};

  EXPECT_TRUE(interpolate(1.0, times, data).isApprox(expected, 1e-10));
}

// 验证只有一个样本点时插值会退化为常值。
TEST(LinearInterpolationTest, InterpolateSingleElement) {
  std::array<double, 1> times = {0.0};
  std::array<Vector<double, 2>, 1> data = {Vector<double, 2>{3.0, 5.0}};

  Vector<double, 2> r = interpolate(0.5, times, data);
  EXPECT_TRUE(r.isApprox(data[0], 1e-10));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
