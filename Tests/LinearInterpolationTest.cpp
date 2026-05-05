/**
 * @file LinearInterpolationTest.cpp
 * @brief 线性插值工具测试：findIndex、findInterval、timeSegment、interpolate。
 */
#include <gtest/gtest.h>
#include "LinearInterpolation.hpp"
#include "Types.hpp"
#include <array>
#include <cmath>

using namespace LinearInterpolation;

TEST(LinearInterpolationTest, FindIndexInTimeArray)
{
    std::array<double, 5> times = {0.0, 1.0, 2.0, 3.0, 4.0};

    EXPECT_EQ(findIndexInTimeArray(times, -0.1), 0u);
    EXPECT_EQ(findIndexInTimeArray(times, 0.0), 0u);
    EXPECT_EQ(findIndexInTimeArray(times, 0.5), 1u);
    EXPECT_EQ(findIndexInTimeArray(times, 1.0), 1u);
    EXPECT_EQ(findIndexInTimeArray(times, 2.5), 3u);
    EXPECT_EQ(findIndexInTimeArray(times, 4.0), 4u);
    EXPECT_EQ(findIndexInTimeArray(times, 5.0), 5u);
}

TEST(LinearInterpolationTest, FindIntervalInTimeArray)
{
    std::array<double, 5> times = {0.0, 1.0, 2.0, 3.0, 4.0};

    EXPECT_EQ(findIntervalInTimeArray(times, -0.1), -1);
    EXPECT_EQ(findIntervalInTimeArray(times, 0.0), -1);
    EXPECT_EQ(findIntervalInTimeArray(times, 0.5), 0);
    EXPECT_EQ(findIntervalInTimeArray(times, 1.0), 0);
    EXPECT_EQ(findIntervalInTimeArray(times, 1.5), 1);
    EXPECT_EQ(findIntervalInTimeArray(times, 3.5), 3);
    EXPECT_EQ(findIntervalInTimeArray(times, 4.0), 3);
    EXPECT_EQ(findIntervalInTimeArray(times, 5.0), 4);
}

TEST(LinearInterpolationTest, TimeSegment)
{
    std::array<double, 4> times = {0.0, 1.0, 2.0, 3.0};

    auto [idx0, alpha0] = timeSegment(0.5, times);
    EXPECT_EQ(idx0, 0);
    EXPECT_DOUBLE_EQ(alpha0, 0.5); // alpha = (t1 - t) / (t1 - t0) = (1 - 0.5)/1 = 0.5

    auto [idx1, alpha1] = timeSegment(1.5, times);
    EXPECT_EQ(idx1, 1);
    EXPECT_DOUBLE_EQ(alpha1, 0.5);

    auto [idx2, alpha2] = timeSegment(0.0, times);
    EXPECT_EQ(idx2, 0);
    EXPECT_DOUBLE_EQ(alpha2, 1.0);

    auto [idx3, alpha3] = timeSegment(3.0, times);
    EXPECT_GE(idx3, 0);
    EXPECT_DOUBLE_EQ(alpha3, 0.0);
}

TEST(LinearInterpolationTest, InterpolateVector)
{
    std::array<double, 3> times = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> data;
    data[0] << 0.0, 0.0;
    data[1] << 1.0, 2.0;
    data[2] << 2.0, 4.0;

    Eigen::Vector2d r0 = interpolate(0.0, times, data);
    EXPECT_TRUE(r0.isApprox(data[0], 1e-10));

    Eigen::Vector2d r1 = interpolate(1.0, times, data);
    EXPECT_TRUE(r1.isApprox(data[1], 1e-10));

    Eigen::Vector2d r05 = interpolate(0.5, times, data);
    Eigen::Vector2d expected05;
    expected05 << 0.5, 1.0;
    EXPECT_TRUE(r05.isApprox(expected05, 1e-10));

    Eigen::Vector2d r15 = interpolate(1.5, times, data);
    Eigen::Vector2d expected15;
    expected15 << 1.5, 3.0;
    EXPECT_TRUE(r15.isApprox(expected15, 1e-10));
}

TEST(LinearInterpolationTest, InterpolateOutOfRangeClampsToBoundary)
{
    std::array<double, 3> times = {0.0, 1.0, 2.0};
    std::array<Eigen::Vector2d, 3> data;
    data[0] << -1.0, 2.0;
    data[1] << 0.0, 4.0;
    data[2] << 3.0, 8.0;

    EXPECT_TRUE(interpolate(-10.0, times, data).isApprox(data[0], 1e-10));
    EXPECT_TRUE(interpolate(10.0, times, data).isApprox(data[2], 1e-10));
}

TEST(LinearInterpolationTest, InterpolateMatrix)
{
    std::array<double, 2> times = {0.0, 2.0};
    std::array<Eigen::Matrix2d, 2> data;
    data[0] << 1.0, 2.0,
               3.0, 4.0;
    data[1] << 5.0, 6.0,
               7.0, 8.0;

    Eigen::Matrix2d expected;
    expected << 3.0, 4.0,
                5.0, 6.0;

    EXPECT_TRUE(interpolate(1.0, times, data).isApprox(expected, 1e-10));
}

TEST(LinearInterpolationTest, InterpolateSingleElement)
{
    std::array<double, 1> times = {0.0};
    std::array<Eigen::Vector2d, 1> data;
    data[0] << 3.0, 5.0;

    Eigen::Vector2d r = interpolate(0.5, times, data);
    EXPECT_TRUE(r.isApprox(data[0], 1e-10));
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
