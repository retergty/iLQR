/**
 * @file RiccatiModification.hpp
 * @brief Riccati 方程修正项：哈密顿量 Hessian、LLT 分解和状态代价修正等。
 */
#pragma once

#include "CholeskyDecomposition.hpp"
#include "Types.hpp"

/**
 * @brief 单节点 Riccati 修正：时间、状态代价修正 deltaQm、哈密顿量 Hessian
 * Hm 及其 Cholesky 分解。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
struct RiccatiModification {
  /** @brief 该节点时间。 */
  Scalar time_ = 0.0;

  /** @brief 状态代价的 Riccati 修正矩阵（如 Hessian 修正等）。 */
  Matrix<Scalar, XDim, XDim> deltaQm_;

  /** @brief 哈密顿量对控制的 Hessian 矩阵 Hm。 */
  Matrix<Scalar, UDim, UDim> hamiltonianHessian_;

  /** @brief Hm 的 Cholesky 分解，用于求解反馈增益和前馈项。 */
  CholeskyDecomposition<Scalar, UDim> HmLLT_;
};
