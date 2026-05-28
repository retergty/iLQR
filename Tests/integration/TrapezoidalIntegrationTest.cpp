/**
 * @file TrapezoidalIntegrationTest.cpp
 * @brief 梯形积分测试：验证标量轨迹积分在简单几何曲线上的结果。
 */
#include <gtest/gtest.h>

#include <array>

#include "Integration/TrapezoidalIntegration.hpp"

TEST(TrapezoidalIntegrationTest, Rectangle) {
  constexpr double width = 8.0;
  constexpr double height = 2.0;
  constexpr size_t numIntervals = 100;
  constexpr double dx = width / static_cast<double>(numIntervals);

  std::array<double, numIntervals + 1> xTrajectory{};
  std::array<double, numIntervals + 1> yTrajectory{};
  yTrajectory.fill(height);

  for (size_t i = 0; i <= numIntervals; ++i) {
    xTrajectory[i] = static_cast<double>(i) * dx;
  }

  const double area = trapezoidalIntegration(xTrajectory, yTrajectory, 0.0);
  EXPECT_NEAR(area, width * height, 1e-12);
}

TEST(TrapezoidalIntegrationTest, RightTriangle) {
  constexpr double width = 8.0;
  constexpr double height = 2.0;
  constexpr size_t numIntervals = 10;
  constexpr double dx = width / static_cast<double>(numIntervals);
  constexpr double dy = height / static_cast<double>(numIntervals);

  std::array<double, numIntervals + 1> xTrajectory{};
  std::array<double, numIntervals + 1> yTrajectory{};

  for (size_t i = 0; i <= numIntervals; ++i) {
    xTrajectory[i] = static_cast<double>(i) * dx;
    yTrajectory[i] = static_cast<double>(i) * dy;
  }

  const double area = trapezoidalIntegration(xTrajectory, yTrajectory, 0.0);
  EXPECT_NEAR(area, width * height / 2.0, 1e-12);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
