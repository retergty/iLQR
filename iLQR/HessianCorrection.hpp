/**
 * @file HessianCorrection.hpp
 * @brief Hessian 修正策略：对角平移等，保证 Riccati 递推数值稳定性。
 */
#pragma once

#include "matrix/Types.hpp"

/**
 * @brief Hessian 矩阵修正策略枚举（当前支持对角平移 DIAGONAL_SHIFT）。
 */
enum class HessianCorrectionStrategy {
  DIAGONAL_SHIFT
  // CHOLESKY_MODIFICATION,
  // EIGENVALUE_MODIFICATION,
  // GERSHGORIN_MODIFICATION
};

/**
 * @brief 按给定策略对 Hessian 矩阵进行修正，使最小特征值不低于 minEigenvalue。
 * @param [in] strategy 修正策略（如对角平移）。
 * @param [in,out] matrix 待修正的 Hessian 矩阵（原地修改）。
 * @param [in] minEigenvalue 修正后期望的最小特征值，默认 1e-6。
 */
template <typename Scalar, int Dimisions>
void shiftHessian(HessianCorrectionStrategy strategy,
                  Matrix<Scalar, Dimisions, Dimisions>& matrix,
                  Scalar minEigenvalue = 1e-6) {
  switch (strategy) {
    case HessianCorrectionStrategy::DIAGONAL_SHIFT: {
      for (int i = 0; i < Dimisions; ++i) {
        matrix(i, i) += minEigenvalue;
      }
      break;
    }
  }
}
