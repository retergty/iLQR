/**
 * @file RiccatiTest.cpp
 * @brief 离散时间 Riccati 方程测试：单步递推的具体数值与对称性。
 */
#include <gtest/gtest.h>

#include "Approximation/LinearApproximation.hpp"
#include "Approximation/QuadraticApproximation.hpp"
#include "ModelData/ModelData.hpp"
#include "RiccatiEquations/DiscreteTimeRiccatiEquations.hpp"
#include "RiccatiEquations/RiccatiModification.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"

// 验证 reduced form 下单步 Riccati 递推与标量参考值一致。
TEST(RiccatiTest, ReducedFormOneStepMatchesScalarReference) {
  using Scalar = double;
  const int XDim = 1;
  const int UDim = 1;

  ModelData<Scalar, XDim, UDim> data;
  data.dynamics.dfdx = {2.0};
  data.dynamics.dfdu = {0.5};
  data.dynamics.f.setZero();
  data.cost.dfdxx = {7.0};
  data.cost.dfdux = {11.0};
  data.cost.dfduu = {13.0};
  data.cost.dfdx = {17.0};
  data.cost.dfdu = {19.0};
  data.cost.f = 23.0;

  RiccatiModification<Scalar, XDim, UDim> mod;
  mod.deltaQm_ = {29.0};

  Matrix<Scalar, XDim, XDim> SmNext{3.0};
  Vector<Scalar, XDim> SvNext{5.0};
  Scalar sNext = 31.0;
  mod.hamiltonianHessian_ = data.cost.dfduu + data.dynamics.dfdu.transpose() *
                                                  SmNext * data.dynamics.dfdu;
  ASSERT_TRUE(mod.HmLLT_.Decomposition(mod.hamiltonianHessian_));

  DiscreteTimeRiccatiEquations<Scalar, XDim, UDim> riccati;

  Matrix<Scalar, UDim, XDim> Km;
  Vector<Scalar, UDim> Lv;
  Matrix<Scalar, XDim, XDim> Sm;
  Vector<Scalar, XDim> Sv;
  Scalar s;

  riccati.computeMap(data, mod, SmNext, SvNext, sNext, Km, Lv, Sm, Sv, s);

  EXPECT_NEAR(Km(0, 0), -56.0 / 55.0, 1e-12);
  EXPECT_NEAR(Lv(0), -86.0 / 55.0, 1e-12);
  EXPECT_NEAR(Sm(0, 0), 1856.0 / 55.0, 1e-12);
  EXPECT_NEAR(Sv(0), 281.0 / 55.0, 1e-12);
  EXPECT_NEAR(s, 4091.0 / 110.0, 1e-12);
}

// 验证零动力学情形下单步 Riccati 递推退化为终端代价形式。
TEST(RiccatiTest, OneStep_ZeroDynamicsTerminalLike) {
  using Scalar = double;
  const int XDim = 2;
  const int UDim = 2;

  ModelData<Scalar, XDim, UDim> data;
  data.dynamics.dfdx.setZero();
  data.dynamics.dfdu = Matrix<Scalar, XDim, UDim>::Identity();
  data.dynamics.f.setZero();
  data.cost.dfdxx = Matrix<Scalar, XDim, XDim>::Identity();
  data.cost.dfdux.setZero();
  data.cost.dfduu = Matrix<Scalar, UDim, UDim>::Identity();
  data.cost.dfdx.setZero();
  data.cost.dfdu.setZero();
  data.cost.f = 0;

  RiccatiModification<Scalar, XDim, UDim> mod;
  mod.deltaQm_.setZero();

  Matrix<Scalar, XDim, XDim> SmNext = Matrix<Scalar, XDim, XDim>::Identity();
  Vector<Scalar, XDim> SvNext = Vector<Scalar, XDim>::Zero();
  Scalar sNext = 0;
  mod.hamiltonianHessian_ = data.cost.dfduu + data.dynamics.dfdu.transpose() *
                                                  SmNext * data.dynamics.dfdu;
  ASSERT_TRUE(mod.HmLLT_.Decomposition(mod.hamiltonianHessian_));

  DiscreteTimeRiccatiEquations<Scalar, XDim, UDim> riccati;

  Matrix<Scalar, UDim, XDim> Km;
  Vector<Scalar, UDim> Lv;
  Matrix<Scalar, XDim, XDim> Sm;
  Vector<Scalar, XDim> Sv;
  Scalar s;

  riccati.computeMap(data, mod, SmNext, SvNext, sNext, Km, Lv, Sm, Sv, s);

  EXPECT_TRUE(Sm.isApprox(Matrix<Scalar, XDim, XDim>::Identity(), 1e-10));
  EXPECT_TRUE(Sv.isApprox(Vector<Scalar, XDim>::Zero(), 1e-10));
  EXPECT_DOUBLE_EQ(s, 0.0);
}

// 验证 non-reduced form 下单步 Riccati 递推与标量参考值一致。
TEST(RiccatiTest, NonReducedFormOneStepMatchesScalarReference) {
  using Scalar = double;
  const int XDim = 1;
  const int UDim = 1;

  ModelData<Scalar, XDim, UDim> data;
  data.dynamics.dfdx = {2.0};
  data.dynamics.dfdu = {0.5};
  data.dynamics.f.setZero();
  data.cost.dfdxx = {7.0};
  data.cost.dfdux = {11.0};
  data.cost.dfduu = {13.0};
  data.cost.dfdx = {17.0};
  data.cost.dfdu = {19.0};
  data.cost.f = 23.0;

  RiccatiModification<Scalar, XDim, UDim> mod;
  mod.deltaQm_ = {29.0};

  Matrix<Scalar, XDim, XDim> SmNext{3.0};
  Vector<Scalar, XDim> SvNext{5.0};
  Scalar sNext = 31.0;
  mod.hamiltonianHessian_ = data.cost.dfduu + data.dynamics.dfdu.transpose() *
                                                  SmNext * data.dynamics.dfdu;
  ASSERT_TRUE(mod.HmLLT_.Decomposition(mod.hamiltonianHessian_));

  DiscreteTimeRiccatiEquations<Scalar, XDim, UDim> riccati;

  Matrix<Scalar, UDim, XDim> Km;
  Vector<Scalar, UDim> Lv;
  Matrix<Scalar, XDim, XDim> Sm;
  Vector<Scalar, XDim> Sv;
  Scalar s;

  riccati.computeMap(data, mod, SmNext, SvNext, sNext, Km, Lv, Sm, Sv, s);

  EXPECT_NEAR(Km(0, 0), -56.0 / 55.0, 1e-12);
  EXPECT_NEAR(Lv(0), -86.0 / 55.0, 1e-12);
  EXPECT_NEAR(Sm(0, 0), 1856.0 / 55.0, 1e-12);
  EXPECT_NEAR(Sv(0), 281.0 / 55.0, 1e-12);
  EXPECT_NEAR(s, 4091.0 / 110.0, 1e-12);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
