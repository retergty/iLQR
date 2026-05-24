//
// 由 rgrandia 于 26.02.20 创建。
//

#pragma once

#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "Cost.hpp"
#include "LinearApproximation.hpp"
#include "LinearStateConstraint.hpp"
#include "LinearStateInputConstraint.hpp"
#include "LinearSystemDynamics.hpp"
#include "QpSolverTypes.hpp"
#include "QpTrajectories.hpp"
#include "QuadraticApproximation.hpp"
#include "QuadraticStateCost.hpp"
#include "StateConstraint.hpp"
#include "StateInputConstraint.hpp"
#include "Types.hpp"

namespace test_tools {

template <typename Scalar, int Dim>
Matrix<Scalar, Dim, Dim> getRandomPositiveDefiniteMatrix() {
  Matrix<Scalar, Dim, Dim> matrix = Matrix<Scalar, Dim, Dim>::Random();
  return matrix.transpose() * matrix + Matrix<Scalar, Dim, Dim>::Identity();
}

template <typename Scalar, int XDim, int UDim>
ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> getRandomCost() {
  static_assert(XDim >= 0, "XDim must be a fixed non-negative dimension.");
  static_assert(UDim >= 0, "UDim must be a fixed non-negative dimension.");

  Matrix<Scalar, XDim + UDim, XDim + UDim> hessian =
      Matrix<Scalar, XDim + UDim, XDim + UDim>::Random();
  hessian = hessian.transpose() * hessian;

  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> cost;
  cost.dfdxx = hessian.template topLeftCorner<XDim, XDim>();
  cost.dfdx = Vector<Scalar, XDim>::Random();
  if constexpr (UDim > 0) {
    cost.dfdux = hessian.template bottomLeftCorner<UDim, XDim>();
    cost.dfduu = hessian.template bottomRightCorner<UDim, UDim>();
    cost.dfdu = Vector<Scalar, UDim>::Random();
  }
  cost.f = std::rand() / static_cast<Scalar>(RAND_MAX);
  return cost;
}

template <typename Scalar, int XDim, int UDim, int ArrayLen>
inline std::unique_ptr<StateInputCost<Scalar, XDim, UDim, ArrayLen>>
getiLQRCost(
    const ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& cost) {
  return std::make_unique<
      QuadraticStateInputCost<Scalar, XDim, UDim, ArrayLen>>(
      cost.dfdxx, cost.dfduu, cost.dfdux);
}

template <typename Scalar, int XDim, int ArrayLen>
inline std::unique_ptr<StateCost<Scalar, XDim, ArrayLen>> getiLQRStateCost(
    const ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>& costFinal) {
  return std::make_unique<QuadraticStateCost<Scalar, XDim, ArrayLen>>(
      costFinal.dfdxx);
}

template <typename Scalar, int XDim, int UDim>
VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
getRandomDynamics() {
  static_assert(XDim >= 0, "XDim must be a fixed non-negative dimension.");
  static_assert(UDim >= 0, "UDim must be a fixed non-negative dimension.");

  VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> dynamics;
  dynamics.dfdx = Matrix<Scalar, XDim, XDim>::Random();
  if constexpr (UDim > 0) {
    dynamics.dfdu = Matrix<Scalar, XDim, UDim>::Random();
  }
  dynamics.f = Vector<Scalar, XDim>::Random();
  return dynamics;
}

template <typename Scalar, int XDim, int UDim>
inline std::unique_ptr<LinearSystemDynamics<Scalar, XDim, UDim>>
getiLQRDynamics(const VectorFunctionLinearApproximation<Scalar, XDim, XDim,
                                                        UDim>& dynamics) {
  return std::make_unique<LinearSystemDynamics<Scalar, XDim, UDim>>(
      dynamics.dfdx, dynamics.dfdu);
}

template <typename Scalar, int XDim, int UDim, int CDim>
inline VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim>
getRandomConstraints() {
  static_assert(XDim >= 0, "XDim must be a fixed non-negative dimension.");
  static_assert(UDim >= 0, "UDim must be a fixed non-negative dimension.");
  static_assert(CDim >= 0, "CDim must be a fixed non-negative dimension.");

  VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim> constraints;
  constraints.dfdx = Matrix<Scalar, CDim, XDim>::Random();
  if constexpr (UDim > 0) {
    constraints.dfdu = Matrix<Scalar, CDim, UDim>::Random();
  }
  constraints.f = Vector<Scalar, CDim>::Random();
  return constraints;
}

template <typename Scalar, int XDim, int UDim, int CDim>
inline std::unique_ptr<StateInputConstraint<Scalar, XDim, UDim, CDim>>
getiLQRConstraints(
    const VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim>&
        stateInputConstraints) {
  return std::make_unique<LinearStateInputConstraint<Scalar, XDim, UDim, CDim>>(
      stateInputConstraints.f, stateInputConstraints.dfdx,
      stateInputConstraints.dfdu);
}

template <typename Scalar, int XDim, int CDim>
inline std::unique_ptr<StateConstraint<Scalar, XDim, CDim>>
getiLQRStateOnlyConstraints(
    const VectorFunctionLinearApproximation<Scalar, CDim, XDim, 0>&
        stateOnlyConstraints) {
  return std::make_unique<LinearStateConstraint<Scalar, XDim, CDim>>(
      stateOnlyConstraints.f, stateOnlyConstraints.dfdx);
}

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
qp_solver::ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>
getRandomTrajectory(Scalar dt) {
  qp_solver::ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> trajectory;
  for (size_t k = 0; k < PredictLength + 1; ++k) {
    trajectory.timeTrajectory[k] = static_cast<Scalar>(k) * dt;
    trajectory.stateTrajectory[k] = Vector<Scalar, XDim>::Random();
  }
  for (size_t k = 0; k < PredictLength; ++k) {
    trajectory.inputTrajectory[k] = Vector<Scalar, UDim>::Random();
  }
  return trajectory;
}

template <typename Scalar, int XDim, int UDim, int CDim>
inline std::vector<qp_solver::LinearQuadraticStage<Scalar, XDim, UDim>>
generateRandomLqProblem(int N) {
  using LinearQuadraticStage_t =
      qp_solver::LinearQuadraticStage<Scalar, XDim, UDim>;

  std::vector<LinearQuadraticStage_t> lqProblem;
  lqProblem.reserve(N + 1);

  for (int k = 0; k < N; ++k) {
    LinearQuadraticStage_t stage;
    const auto constraints = getRandomConstraints<Scalar, XDim, UDim, CDim>();
    stage.cost = getRandomCost<Scalar, XDim, UDim>();
    stage.dynamics = getRandomDynamics<Scalar, XDim, UDim>();
    stage.constraints.f = constraints.f;
    stage.constraints.dfdx = constraints.dfdx;
    if constexpr (UDim > 0) {
      stage.constraints.dfdu = constraints.dfdu;
    }
    lqProblem.emplace_back(std::move(stage));
  }

  LinearQuadraticStage_t terminalStage;
  const auto terminalCost = getRandomCost<Scalar, XDim, 0>();
  const auto terminalConstraints =
      getRandomConstraints<Scalar, XDim, UDim, CDim>();
  terminalStage.cost.setZero();
  terminalStage.cost.dfdxx = terminalCost.dfdxx;
  terminalStage.cost.dfdx = terminalCost.dfdx;
  terminalStage.cost.f = terminalCost.f;
  terminalStage.dynamics.setZero();
  terminalStage.constraints.f = terminalConstraints.f;
  terminalStage.constraints.dfdx = terminalConstraints.dfdx;
  if constexpr (UDim > 0) {
    terminalStage.constraints.dfdu = terminalConstraints.dfdu;
  }
  lqProblem.emplace_back(std::move(terminalStage));

  return lqProblem;
}

/** 检查 QP 可行性和数值条件。 */
template <typename Scalar, int DecisionDim, int ConstraintDim>
inline bool isQpFeasible(
    const ScalarFunctionQuadraticApproximation<Scalar, DecisionDim, 0>& qpCost,
    const VectorFunctionLinearApproximation<Scalar, ConstraintDim, DecisionDim,
                                            0>& qpConstraints) {
  const auto& H = qpCost.dfdxx;
  const auto& A = qpConstraints.dfdx;

  // 代价必须凸。
  Eigen::LDLT<Matrix<Scalar, DecisionDim, DecisionDim>> ldlt(H);
  if (!(ldlt.vectorD().array() > Scalar(0.0)).all()) {
    std::cerr << "H is not positive definite\n";
    return false;
  }

  if (A.rows() == 0) {
    return true;
  }

  // 约束可行性。
  Eigen::JacobiSVD<Matrix<Scalar, ConstraintDim, DecisionDim>> svd(A);
  const auto conditionNumber =
      svd.singularValues()(0) / svd.singularValues().tail(1)(0);
  if (svd.rank() != A.rows()) {
    std::cerr << "A is not full row-rank\n";
    return false;
  } else if (conditionNumber > 1e6) {
    std::cerr << "A is ill-conditioned\n";
    return false;
  }

  return true;
}

}  // namespace test_tools
