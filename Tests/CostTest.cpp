/**
 * @file CostTest.cpp
 * @brief 代价模块测试：匹配当前项目中的 QuadraticStateCost /
 * QuadraticStateInputCost 接口。
 */
#include <gtest/gtest.h>

#include <array>

#include "Cost/QuadraticStateCost.hpp"
#include "Types.hpp"

class QuadraticCostTest : public testing::Test {
 protected:
  using Scalar = double;
  static constexpr int XDim = 2;
  static constexpr int UDim = 1;
  static constexpr int ArrayLength = 1;

  using StateMatrix = Matrix<Scalar, XDim, XDim>;
  using InputMatrix = Matrix<Scalar, UDim, UDim>;
  using InputStateMatrix = Matrix<Scalar, UDim, XDim>;
  using StateVector = Vector<Scalar, XDim>;
  using InputVector = Vector<Scalar, UDim>;

  static constexpr Scalar kPrecision = 1e-12;

  QuadraticCostTest() {
    Q_ = {{2.0, 1.0}, {1.0, 2.0}};
    Qf_ = Matrix<Scalar, XDim, XDim>::Identity();
    R_ = {2.0};
    P_ = {1.0, 1.0};

    t_ = 0.0;

    xNominal_ = {1.0, -2.0};
    uNominal_ = {0.5};
    x_ = {3.0, -1.0};
    u_ = {2.0};

    timeTrajectory_ = {t_};
    stateTrajectory_ = {xNominal_};
    inputTrajectory_ = {uNominal_};

    const StateVector dx = x_ - xNominal_;
    const InputVector du = u_ - uNominal_;
    expectedCost_ =
        0.5 * dx.dot(Q_ * dx) + 0.5 * du.dot(R_ * du) + du.dot(P_ * dx);
    expectedFinalCost_ = 0.5 * dx.dot(Qf_ * dx);
  }

  StateMatrix Q_;
  StateMatrix Qf_;
  InputMatrix R_;
  InputStateMatrix P_;
  Scalar t_{0.0};
  StateVector x_;
  InputVector u_;
  StateVector xNominal_;
  InputVector uNominal_;
  std::array<Scalar, ArrayLength> timeTrajectory_;
  std::array<StateVector, ArrayLength> stateTrajectory_;
  std::array<InputVector, ArrayLength> inputTrajectory_;
  Scalar expectedCost_{0.0};
  Scalar expectedFinalCost_{0.0};
};

TEST_F(QuadraticCostTest, StateInputCostValue) {
  QuadraticStateInputCost<Scalar, XDim, UDim, ArrayLength> costFunction(Q_, R_,
                                                                        P_, 0);

  const auto cost = costFunction.getValue(t_, x_, u_, timeTrajectory_,
                                          stateTrajectory_, inputTrajectory_);
  EXPECT_NEAR(cost, expectedCost_, kPrecision);
}

TEST_F(QuadraticCostTest, QuadraticCostsStoreExplicitCostNumber) {
  QuadraticStateInputCost<Scalar, XDim, UDim, ArrayLength> stateInputCost(
      Q_, R_, P_, 17);
  QuadraticStateCost<Scalar, XDim, ArrayLength> stateCost(Qf_, 23);

  EXPECT_EQ(stateInputCost.number, 17);
  EXPECT_EQ(stateCost.number, 23);
}

TEST_F(QuadraticCostTest, StateInputCostApproximation) {
  QuadraticStateInputCost<Scalar, XDim, UDim, ArrayLength> costFunction(Q_, R_,
                                                                        P_, 0);

  const auto approximation = costFunction.getQuadraticApproximation(
      t_, x_, u_, timeTrajectory_, stateTrajectory_, inputTrajectory_);

  const StateVector dx = x_ - xNominal_;
  const InputVector du = u_ - uNominal_;

  EXPECT_NEAR(approximation.f, expectedCost_, kPrecision);
  EXPECT_TRUE(
      approximation.dfdx.isApprox(Q_ * dx + P_.transpose() * du, kPrecision));
  EXPECT_TRUE(approximation.dfdu.isApprox(R_ * du + P_ * dx, kPrecision));
  EXPECT_TRUE(approximation.dfdxx.isApprox(Q_, kPrecision));
  EXPECT_TRUE(approximation.dfdux.isApprox(P_, kPrecision));
  EXPECT_TRUE(approximation.dfduu.isApprox(R_, kPrecision));
}

TEST_F(QuadraticCostTest, StateCostValue) {
  QuadraticStateCost<Scalar, XDim, ArrayLength> costFunction(Qf_, 0);

  const auto cost =
      costFunction.getValue(t_, x_, timeTrajectory_, stateTrajectory_);
  EXPECT_NEAR(cost, expectedFinalCost_, kPrecision);
}

TEST_F(QuadraticCostTest, StateCostApproximation) {
  QuadraticStateCost<Scalar, XDim, ArrayLength> costFunction(Qf_, 0);

  const auto approximation = costFunction.getQuadraticApproximation(
      t_, x_, timeTrajectory_, stateTrajectory_);

  const StateVector dx = x_ - xNominal_;

  EXPECT_NEAR(approximation.f, expectedFinalCost_, kPrecision);
  EXPECT_TRUE(approximation.dfdx.isApprox(Qf_ * dx, kPrecision));
  EXPECT_TRUE(approximation.dfdxx.isApprox(Qf_, kPrecision));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}