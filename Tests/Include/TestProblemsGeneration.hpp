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

#include "Approximation/LinearApproximation.hpp"
#include "Approximation/QuadraticApproximation.hpp"
#include "Constraint/LinearStateConstraint.hpp"
#include "Constraint/LinearStateInputConstraint.hpp"
#include "Constraint/StateConstraint.hpp"
#include "Constraint/StateInputConstraint.hpp"
#include "Cost/Cost.hpp"
#include "Cost/QuadraticStateCost.hpp"
#include "Dynamics/LinearSystemDynamics.hpp"
#include "Tests/Include/MatrixEigenConversion.hpp"
#include "Tests/Include/QpSolverTypes.hpp"
#include "Tests/Include/QpTrajectories.hpp"

namespace test_tools {

using matrix_eigen_conversion::fromEigenMatrix;
using matrix_eigen_conversion::fromEigenVector;
using matrix_eigen_conversion::toEigenMatrix;
using matrix_eigen_conversion::toEigenVector;

template <typename Scalar, int Dim>
Eigen::Matrix<Scalar, Dim, Dim> getRandomPositiveDefiniteMatrix() {
  Eigen::Matrix<Scalar, Dim, Dim> matrix =
      Eigen::Matrix<Scalar, Dim, Dim>::Random();
  return matrix.transpose() * matrix +
         Eigen::Matrix<Scalar, Dim, Dim>::Identity();
}

template <typename Scalar, int XDim, int UDim>
ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> getRandomCost() {
  static_assert(XDim >= 0, "XDim must be a fixed non-negative dimension.");
  static_assert(UDim >= 0, "UDim must be a fixed non-negative dimension.");

  Eigen::Matrix<Scalar, XDim + UDim, XDim + UDim> hessian =
      Eigen::Matrix<Scalar, XDim + UDim, XDim + UDim>::Random();
  hessian = hessian.transpose() * hessian;

  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> cost;
  cost.dfdxx = fromEigenMatrix(hessian.template topLeftCorner<XDim, XDim>());
  cost.dfdx = fromEigenVector(Eigen::Vector<Scalar, XDim>::Random());
  if constexpr (UDim > 0) {
    cost.dfdux =
        fromEigenMatrix(hessian.template bottomLeftCorner<UDim, XDim>());
    cost.dfduu =
        fromEigenMatrix(hessian.template bottomRightCorner<UDim, UDim>());
    cost.dfdu = fromEigenVector(Eigen::Vector<Scalar, UDim>::Random());
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
      cost.dfdxx, cost.dfduu, cost.dfdux, 0);
}

template <typename Scalar, int XDim, int UDim, int ArrayLen>
inline std::unique_ptr<StateCost<Scalar, XDim, UDim, ArrayLen>>
getiLQRStateCost(
    const ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>& costFinal) {
  return std::make_unique<QuadraticStateCost<Scalar, XDim, UDim, ArrayLen>>(
      costFinal.dfdxx, 0);
}

template <typename Scalar, int XDim, int UDim>
VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
getRandomDynamics() {
  static_assert(XDim >= 0, "XDim must be a fixed non-negative dimension.");
  static_assert(UDim >= 0, "UDim must be a fixed non-negative dimension.");

  VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> dynamics;
  dynamics.dfdx = fromEigenMatrix(Eigen::Matrix<Scalar, XDim, XDim>::Random());
  if constexpr (UDim > 0) {
    dynamics.dfdu =
        fromEigenMatrix(Eigen::Matrix<Scalar, XDim, UDim>::Random());
  }
  dynamics.f = fromEigenVector(Eigen::Vector<Scalar, XDim>::Random());
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
  constraints.dfdx =
      fromEigenMatrix(Eigen::Matrix<Scalar, CDim, XDim>::Random());
  if constexpr (UDim > 0) {
    constraints.dfdu =
        fromEigenMatrix(Eigen::Matrix<Scalar, CDim, UDim>::Random());
  }
  constraints.f = fromEigenVector(Eigen::Vector<Scalar, CDim>::Random());
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
    trajectory.stateTrajectory[k] = Eigen::Vector<Scalar, XDim>::Random();
  }
  for (size_t k = 0; k < PredictLength; ++k) {
    trajectory.inputTrajectory[k] = Eigen::Vector<Scalar, UDim>::Random();
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
    stage.constraints.f = toEigenVector(constraints.f);
    stage.constraints.dfdx = toEigenMatrix(constraints.dfdx);
    if constexpr (UDim > 0) {
      stage.constraints.dfdu = toEigenMatrix(constraints.dfdu);
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
  terminalStage.constraints.f = toEigenVector(terminalConstraints.f);
  terminalStage.constraints.dfdx = toEigenMatrix(terminalConstraints.dfdx);
  if constexpr (UDim > 0) {
    terminalStage.constraints.dfdu = toEigenMatrix(terminalConstraints.dfdu);
  }
  lqProblem.emplace_back(std::move(terminalStage));

  return lqProblem;
}

/** 检查动态 QP 可行性和数值条件。 */
template <typename Scalar>
inline bool isQpFeasible(
    const qp_solver::QpCostApproximation<Scalar>& qpCost,
    const qp_solver::QpDenseConstraintApproximation<Scalar>& qpConstraints) {
  const auto& H = qpCost.dfdxx;
  const auto& A = qpConstraints.dfdx;

  Eigen::LDLT<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> ldlt(H);
  if (!(ldlt.vectorD().array() > Scalar(0.0)).all()) {
    std::cerr << "H is not positive definite\n";
    return false;
  }

  if (A.rows() == 0) {
    return true;
  }

  Eigen::JacobiSVD<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> svd(
      A);
  const auto conditionNumber =
      svd.singularValues()(0) / svd.singularValues().tail(1)(0);
  if (svd.rank() != A.rows()) {
    std::cerr << "A is not full row-rank\n";
    return false;
  } else if (conditionNumber > Scalar(1e6)) {
    std::cerr << "A is ill-conditioned\n";
    return false;
  }

  return true;
}

}  // namespace test_tools
