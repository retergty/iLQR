#pragma once

#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"

namespace qp_solver {

/** 定义给定阶段的二次代价和线性动力学。 */
template <typename Scalar, int XDim, int UDim>
struct LinearQuadraticStage {
  using CostApproximation_t =
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>;
  using DynamicsApproximation_t =
      VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>;
  using ConstraintApproximation_t =
      VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, XDim, UDim>;
  /** 代价的二次近似。 */
  CostApproximation_t cost;
  /** 动力学的线性近似。 */
  DynamicsApproximation_t dynamics;
  /** 约束的线性近似。 */
  ConstraintApproximation_t constraints;

  LinearQuadraticStage() {
    constraints.f.resize(0);
    constraints.dfdx.resize(0, XDim);
    constraints.dfdu.resize(0, UDim);
  }

  LinearQuadraticStage(CostApproximation_t c, DynamicsApproximation_t d,
                       ConstraintApproximation_t g)
      : cost(std::move(c)), dynamics(std::move(d)), constraints(std::move(g)) {}
};

}  // namespace qp_solver
