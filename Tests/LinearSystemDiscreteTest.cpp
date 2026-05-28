#include <gtest/gtest.h>

#include "Dynamics/LinearSystemDynamics.hpp"
#include "Integration/SensitivityIntegrator.hpp"
#include "math.h"

// 验证 EK2 离散化结果与显式 RK2 公式一致。
TEST(LinearSystemTest, DiscreteEK2Test) {
  Matrix<double, 2, 2> A{{0.0, 1.0}, {0.0, 0.0}};

  Matrix<double, 2, 1> B{0.0, 1.0};

  LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
  EK2DynamicsDiscretizer<double, 2, 1> ek2Discretizer;

  Vector<double, 2> x{1.0, 2.0};

  Vector<double, 1> u{3.0};

  const double t = 0.0;
  const double dt = 0.5;
  const Vector<double, 2> k1 = A * x + B * u;
  const Vector<double, 2> k2 = A * (x + dt * k1) + B * u;
  const Vector<double, 2> expected = x + 0.5 * dt * (k1 + k2);

  const Vector<double, 2> actual =
      ek2Discretizer.discretize(lin_sys, t, x, u, dt);

  EXPECT_TRUE(actual.isApprox(expected, 1e-12));
}

// 验证前向欧拉离散化结果与线性系统的显式离散更新一致。
TEST(LinearSystemTest, DiscreteEulerTest) {
  Matrix<double, 2, 2> A{{1.0, 2.0}, {0.0, -1.0}};

  Matrix<double, 2, 1> B{0.5, 2.0};

  LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
  EulerDynamicsDiscretizer<double, 2, 1> eulerDiscretizer;

  Vector<double, 2> x{1.0, -2.0};

  Vector<double, 1> u{3.0};

  const double dt = 0.25;
  const Vector<double, 2> expected = x + dt * (A * x + B * u);
  const Vector<double, 2> actual =
      eulerDiscretizer.discretize(lin_sys, 1.0, x, u, dt);

  EXPECT_TRUE(actual.isApprox(expected, 1e-12));
}

// 验证 EK4 离散化结果与经典 RK4 公式一致。
TEST(LinearSystemTest, DiscreteEK4Test) {
  Matrix<double, 2, 2> A{{0.0, 1.0}, {-2.0, -3.0}};

  Matrix<double, 2, 1> B{0.0, 1.0};

  LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
  EK4DynamicsDiscretizer<double, 2, 1> ek4Discretizer;

  Vector<double, 2> x{1.0, -1.0};

  Vector<double, 1> u{2.0};

  const double t = 0.0;
  const double dt = 0.1;
  const Vector<double, 2> k1 = A * x + B * u;
  const Vector<double, 2> k2 = A * (x + 0.5 * dt * k1) + B * u;
  const Vector<double, 2> k3 = A * (x + 0.5 * dt * k2) + B * u;
  const Vector<double, 2> k4 = A * (x + dt * k3) + B * u;
  const Vector<double, 2> expected =
      x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

  const Vector<double, 2> actual =
      ek4Discretizer.discretize(lin_sys, t, x, u, dt);

  EXPECT_TRUE(actual.isApprox(expected, 1e-12));
}

// 验证欧拉离散化的一阶敏感度与线性系统解析 Jacobian 一致。
TEST(LinearSystemTest, EulerSensitivityMatchesLinearSystemAnalyticResult) {
  Matrix<double, 2, 2> A{{1.0, 2.0}, {0.0, -1.0}};

  Matrix<double, 2, 1> B{0.5, 2.0};

  LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
  EulerDynamicsDiscretizer<double, 2, 1> eulerDiscretizer;

  Vector<double, 2> x{1.0, -2.0};

  Vector<double, 1> u{3.0};

  const double dt = 0.25;
  const auto approx =
      eulerDiscretizer.sensitivityDiscretize(lin_sys, 1.0, x, u, dt);

  EXPECT_TRUE(
      approx.dfdx.isApprox(Matrix<double, 2, 2>::Identity() + dt * A, 1e-12));
  EXPECT_TRUE(approx.dfdu.isApprox(dt * B, 1e-12));
  EXPECT_TRUE(approx.f.isApprox(x + dt * (A * x + B * u), 1e-12));
}

// 验证 EK2 离散化的一阶敏感度与线性系统二阶展开一致。
TEST(LinearSystemTest, EK2SensitivityMatchesLinearSystemAnalyticResult) {
  Matrix<double, 2, 2> A{{0.0, 1.0}, {-2.0, -3.0}};

  Matrix<double, 2, 1> B{0.0, 1.0};

  LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
  EK2DynamicsDiscretizer<double, 2, 1> ek2Discretizer;

  Vector<double, 2> x{1.0, -1.0};

  Vector<double, 1> u{2.0};

  const double dt = 0.1;
  const auto approx =
      ek2Discretizer.sensitivityDiscretize(lin_sys, 0.0, x, u, dt);
  const Matrix<double, 2, 2> expectedDfdx =
      Matrix<double, 2, 2>::Identity() + dt * A + 0.5 * dt * dt * A * A;
  const Matrix<double, 2, 1> expectedDfdu = dt * B + 0.5 * dt * dt * A * B;
  const Vector<double, 2> expectedF =
      ek2Discretizer.discretize(lin_sys, 0.0, x, u, dt);

  EXPECT_TRUE(approx.dfdx.isApprox(expectedDfdx, 1e-12));
  EXPECT_TRUE(approx.dfdu.isApprox(expectedDfdu, 1e-12));
  EXPECT_TRUE(approx.f.isApprox(expectedF, 1e-12));
}

// 验证 EK4 离散化的一阶敏感度与线性系统四阶展开一致。
TEST(LinearSystemTest, EK4SensitivityMatchesLinearSystemAnalyticResult) {
  Matrix<double, 2, 2> A{{0.0, 1.0}, {-2.0, -3.0}};

  Matrix<double, 2, 1> B{0.0, 1.0};

  LinearSystemDynamics<double, 2, 1> lin_sys(A, B);
  EK4DynamicsDiscretizer<double, 2, 1> ek4Discretizer;

  Vector<double, 2> x{1.0, -1.0};

  Vector<double, 1> u{2.0};

  const double dt = 0.1;
  const auto approx =
      ek4Discretizer.sensitivityDiscretize(lin_sys, 0.0, x, u, dt);
  const Matrix<double, 2, 2> expectedDfdx =
      Matrix<double, 2, 2>::Identity() + dt * A + 0.5 * dt * dt * A * A +
      (dt * dt * dt / 6.0) * A * A * A +
      (dt * dt * dt * dt / 24.0) * A * A * A * A;
  const Matrix<double, 2, 1> expectedDfdu =
      dt * B + 0.5 * dt * dt * A * B + (dt * dt * dt / 6.0) * A * A * B +
      (dt * dt * dt * dt / 24.0) * A * A * A * B;
  const Vector<double, 2> expectedF =
      ek4Discretizer.discretize(lin_sys, 0.0, x, u, dt);

  EXPECT_TRUE(approx.dfdx.isApprox(expectedDfdx, 1e-12));
  EXPECT_TRUE(approx.dfdu.isApprox(expectedDfdu, 1e-12));
  EXPECT_TRUE(approx.f.isApprox(expectedF, 1e-12));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}