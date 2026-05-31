/**
 * @file ModelData.hpp
 * @brief 单时刻最优控制问题的线性化/二次近似数据：动力学与代价的系数。
 */
#pragma once

#include "Approximation/LinearApproximation.hpp"
#include "Approximation/QuadraticApproximation.hpp"

/**
 * @brief 单节点 LQ 模型数据：时间、动力学线性近似（dfdx, dfdu,
 * f）、代价二次近似（dfdxx, dfdux, dfduu, dfdx, dfdu, f）。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
struct ModelData {
  Scalar time = Scalar(0.0);

  /** @brief 动力学线性近似：x_{k+1} ≈ dfdx*dx + dfdu*du + f。 */
  VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> dynamics;

  /** @brief 代价二次近似（标量函数对 (x,u) 的系数）。 */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> cost;
};