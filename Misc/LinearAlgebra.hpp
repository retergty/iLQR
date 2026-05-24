/**
 * @file LinearAlgebra.hpp
 * @brief 线性代数工具：对称正定矩阵的逆的 UUT 分解（用于投影等）。
 */
#pragma once

#include <Eigen/Cholesky>
#include <Eigen/LU>

#include "Types.hpp"

namespace LinearAlgebra {
/**
 * @brief 计算 inv(Am) 的 U*U^T 分解中的上三角因子 U（Am =
 * L*L^T，inv(Am)=inv(L^T)*inv(L)）。
 * @param [in] Am 对称正定方阵。
 * @param [out] AmInvUmUmT 上三角矩阵，满足 AmInvUmUmT * AmInvUmUmT^T =
 * inv(Am)。
 */
template <typename Scalar, int Dimisions>
void computeInverseMatrixUUT(const Matrix<Scalar, Dimisions, Dimisions>& Am,
                             Matrix<Scalar, Dimisions, Dimisions>& AmInvUmUmT) {
  // Am = Lm Lm^T --> inv(Am) = inv(Lm^T) inv(Lm)，其中 Lm^T 为上三角矩阵。
  Eigen::LLT<Matrix<Scalar, Dimisions, Dimisions>> lltOfA(Am);
  AmInvUmUmT.setIdentity(Am.rows(), Am.cols());  // 用于动态尺寸矩阵。
  lltOfA.matrixU().solveInPlace(AmInvUmUmT);
}
}  // namespace LinearAlgebra