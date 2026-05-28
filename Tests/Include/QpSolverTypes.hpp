#pragma once

#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"

namespace qp_solver {

/** QP 测试工具专用的动态二次代价近似。 */
template <typename Scalar>
struct QpCostApproximation {
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dfdxx;
  Eigen::Vector<Scalar, Eigen::Dynamic> dfdx;
  Scalar f{0};
};

/** QP 测试工具专用的完整动态线性约束近似。 */
template <typename Scalar>
struct QpDenseConstraintApproximation {
  Eigen::Vector<Scalar, Eigen::Dynamic> f;
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dfdx;
};

/** QP 测试工具专用的动态线性约束近似。 */
template <typename Scalar, int XDim, int UDim>
struct QpConstraintApproximation {
  Eigen::Vector<Scalar, Eigen::Dynamic> f;
  Eigen::Matrix<Scalar, Eigen::Dynamic, XDim> dfdx;
  Eigen::Matrix<Scalar, Eigen::Dynamic, UDim> dfdu;

  QpConstraintApproximation() {
    f.resize(0);
    dfdx.resize(0, XDim);
    dfdu.resize(0, UDim);
  }
};

/** 定义给定阶段的二次代价和线性动力学。 */
template <typename Scalar, int XDim, int UDim>
struct LinearQuadraticStage {
  using CostApproximation_t =
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>;
  using DynamicsApproximation_t =
      VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>;
  using ConstraintApproximation_t =
      QpConstraintApproximation<Scalar, XDim, UDim>;
  /** 代价的二次近似。 */
  CostApproximation_t cost;
  /** 动力学的线性近似。 */
  DynamicsApproximation_t dynamics;
  /** 约束的线性近似。 */
  ConstraintApproximation_t constraints;

  LinearQuadraticStage() = default;

  LinearQuadraticStage(CostApproximation_t c, DynamicsApproximation_t d,
                       ConstraintApproximation_t g)
      : cost(std::move(c)), dynamics(std::move(d)), constraints(std::move(g)) {}
};

}  // namespace qp_solver
