#include <gtest/gtest.h>

#include <array>

#include "Integration/Observer.hpp"
#include "Integration/TrapezoidalIntegration.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"

// 验证非均匀时间网格上的标量梯形积分结果正确。
TEST(IntegrationUtilitiesTest, TrapezoidalIntegrationScalarNonUniformGrid) {
  std::array<double, 3> times = {0.0, 1.0, 3.0};
  std::array<double, 3> values = {0.0, 2.0, 4.0};

  const double result = trapezoidalIntegration(times, values, 1.0);

  EXPECT_DOUBLE_EQ(result, 8.0);
}

// 验证向量梯形积分按分量独立计算。
TEST(IntegrationUtilitiesTest, TrapezoidalIntegrationVector) {
  std::array<double, 2> times = {0.0, 2.0};
  std::array<Vector<double, 2>, 2> values = {Vector<double, 2>{1.0, -1.0},
                                             Vector<double, 2>{3.0, 5.0}};

  const Vector<double, 2> initial{10.0, 20.0};
  const Vector<double, 2> expected{14.0, 24.0};

  const Vector<double, 2> result =
      trapezoidalIntegration(times, values, initial);

  EXPECT_TRUE(result.isApprox(expected, 1e-10));
}

// 验证 Observer 达到容量后停止写入，并且可通过 clear() 重置。
TEST(IntegrationUtilitiesTest, ObserverStoresUpToCapacityAndCanClear) {
  std::array<double, 2> times = {};
  std::array<Vector<double, 2>, 2> states = {Vector<double, 2>::Zero(),
                                             Vector<double, 2>::Zero()};

  Observer<double, 2> observer(2, states.data(), times.data());

  const Vector<double, 2> x0{1.0, 2.0};
  const Vector<double, 2> x1{3.0, 4.0};
  const Vector<double, 2> x2{5.0, 6.0};

  observer.observe(x0, 0.0);
  observer.observe(x1, 1.0);
  observer.observe(x2, 2.0);

  EXPECT_EQ(observer.getCount(), 2);
  EXPECT_DOUBLE_EQ(times[0], 0.0);
  EXPECT_DOUBLE_EQ(times[1], 1.0);
  EXPECT_TRUE(states[0].isApprox(x0, 1e-10));
  EXPECT_TRUE(states[1].isApprox(x1, 1e-10));

  observer.clear();
  observer.observe(x2, 3.0);

  EXPECT_EQ(observer.getCount(), 1);
  EXPECT_DOUBLE_EQ(times[0], 3.0);
  EXPECT_TRUE(states[0].isApprox(x2, 1e-10));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
