/**
 * @file RiccatiTest.cpp
 * @brief 离散时间 Riccati 方程测试：单步递推的具体数值与对称性。
 */
#include "DiscreteTimeRiccatiEquations.hpp"
#include "LinearApproximation.hpp"
#include "ModelData.hpp"
#include "QuadraticApproximation.hpp"
#include "RiccatiModification.hpp"
#include "Types.hpp"
#include <gtest/gtest.h>

// 验证 reduced form 下单步 Riccati 递推与标量参考值一致。
TEST(RiccatiTest, ReducedFormOneStepMatchesScalarReference) {
  using Scalar = double;
  const int XDim = 1;
  const int UDim = 1;

  ModelData<Scalar, XDim, UDim> data;
  data.dynamics.dfdx << 2.0;
  data.dynamics.dfdu << 0.5;
  data.dynamics.f.setZero();
  data.cost.dfdxx << 7.0;
  data.cost.dfdux << 11.0;
  data.cost.dfduu << 13.0;
  data.cost.dfdx << 17.0;
  data.cost.dfdu << 19.0;
  data.cost.f = 23.0;

  RiccatiModification<Scalar, XDim, UDim> mod;
  mod.deltaQm_ << 29.0;

  Eigen::Matrix<double, 1, 1> SmNext;
  SmNext << 3.0;
  Eigen::Matrix<double, 1, 1> SvNext;
  SvNext << 5.0;
  Scalar sNext = 31.0;

  DiscreteTimeRiccatiEquations<Scalar, XDim, UDim> riccati(true);

  Eigen::Matrix<double, 1, 1> Km;
  Eigen::Matrix<double, 1, 1> Lv;
  Eigen::Matrix<double, 1, 1> Sm;
  Eigen::Matrix<double, 1, 1> Sv;
  Scalar s;

  riccati.computeMap(data, mod, SmNext, SvNext, sNext, Km, Lv, Sm, Sv, s);

  EXPECT_NEAR(Km(0, 0), -14.0, 1e-12);
  EXPECT_NEAR(Lv(0), -21.5, 1e-12);
  EXPECT_NEAR(Sm(0, 0), -148.0, 1e-12);
  EXPECT_NEAR(Sv(0), -274.0, 1e-12);
  EXPECT_NEAR(s, -177.125, 1e-12);
}

// 验证零动力学情形下单步 Riccati 递推退化为终端代价形式。
TEST(RiccatiTest, OneStep_ZeroDynamicsTerminalLike) {
  using Scalar = double;
  const int XDim = 2;
  const int UDim = 2;

  ModelData<Scalar, XDim, UDim> data;
  data.dynamics.dfdx.setZero();
  data.dynamics.dfdu = Eigen::Matrix2d::Identity();
  data.dynamics.f.setZero();
  data.cost.dfdxx = Eigen::Matrix2d::Identity();
  data.cost.dfdux.setZero();
  data.cost.dfduu = Eigen::Matrix2d::Identity();
  data.cost.dfdx.setZero();
  data.cost.dfdu.setZero();
  data.cost.f = 0;

  RiccatiModification<Scalar, XDim, UDim> mod;
  mod.deltaQm_.setZero();

  Eigen::Matrix2d SmNext = Eigen::Matrix2d::Identity();
  Eigen::Vector2d SvNext = Eigen::Vector2d::Zero();
  Scalar sNext = 0;

  DiscreteTimeRiccatiEquations<Scalar, XDim, UDim> riccati(true);

  Eigen::Matrix2d Km;
  Eigen::Vector2d Lv;
  Eigen::Matrix2d Sm;
  Eigen::Vector2d Sv;
  Scalar s;

  riccati.computeMap(data, mod, SmNext, SvNext, sNext, Km, Lv, Sm, Sv, s);

  EXPECT_TRUE(Sm.isApprox(Eigen::Matrix2d::Identity(), 1e-10));
  EXPECT_TRUE(Sv.isApprox(Eigen::Vector2d::Zero(), 1e-10));
  EXPECT_DOUBLE_EQ(s, 0.0);
}

// 验证 non-reduced form 下单步 Riccati 递推与标量参考值一致。
TEST(RiccatiTest, NonReducedFormOneStepMatchesScalarReference) {
  using Scalar = double;
  const int XDim = 1;
  const int UDim = 1;

  ModelData<Scalar, XDim, UDim> data;
  data.dynamics.dfdx << 2.0;
  data.dynamics.dfdu << 0.5;
  data.dynamics.f.setZero();
  data.cost.dfdxx << 7.0;
  data.cost.dfdux << 11.0;
  data.cost.dfduu << 13.0;
  data.cost.dfdx << 17.0;
  data.cost.dfdu << 19.0;
  data.cost.f = 23.0;

  RiccatiModification<Scalar, XDim, UDim> mod;
  mod.deltaQm_ << 29.0;

  Eigen::Matrix<double, 1, 1> SmNext;
  SmNext << 3.0;
  Eigen::Matrix<double, 1, 1> SvNext;
  SvNext << 5.0;
  Scalar sNext = 31.0;

  DiscreteTimeRiccatiEquations<Scalar, XDim, UDim> riccati(false);

  Eigen::Matrix<double, 1, 1> Km;
  Eigen::Matrix<double, 1, 1> Lv;
  Eigen::Matrix<double, 1, 1> Sm;
  Eigen::Matrix<double, 1, 1> Sv;
  Scalar s;

  riccati.computeMap(data, mod, SmNext, SvNext, sNext, Km, Lv, Sm, Sv, s);

  EXPECT_NEAR(Km(0, 0), -14.0, 1e-12);
  EXPECT_NEAR(Lv(0), -21.5, 1e-12);
  EXPECT_NEAR(Sm(0, 0), 2351.0, 1e-12);
  EXPECT_NEAR(Sv(0), 3563.75, 1e-12);
  EXPECT_NEAR(s, 2769.71875, 1e-12);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
