/**
 * @file ChangeOfInputVariables.hpp
 * @brief 输入变量变换：对二次近似做 δu = Pu * δũ，将输入维度由 m 变为 p（Pu 为
 * m×p 矩阵）。
 */
#pragma once

#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "Types.hpp"

/**
 * @brief 对二次近似施加输入变换 δu = Pu * δũ，变换后状态维 n、输入维 p；Pu
 * 为空矩阵时表示零输入。
 * @param [in,out] quadraticApproximation 待变换的二次近似（原地修改）。
 * @param [in] Pu 定义 δũ 范围的矩阵（m×p）。
 */
template <typename Scalar, int XDim, int UDim>
void changeOfInputVariables(ScalarFunctionQuadraticApproximation<
                                Scalar, XDim, UDim>& quadraticApproximation,
                            const Matrix<Scalar, UDim, UDim>& Pu) {
  // P = Pu'*P
  quadraticApproximation.dfdux = Pu.transpose() * quadraticApproximation.dfdux;

  // R = Pu' * R * Pu
  quadraticApproximation.dfduu =
      Pu.transpose() * quadraticApproximation.dfduu * Pu;

  // r = Pu' * r
  quadraticApproximation.dfdu = Pu.transpose() * quadraticApproximation.dfdu;
}

/** 将输入变量变换应用到线性系统。 */
template <typename Scalar, int FDimisions, int XDim, int UDim>
void changeOfInputVariables(
    VectorFunctionLinearApproximation<Scalar, FDimisions, XDim, UDim>&
        linearApproximation,
    const Matrix<Scalar, UDim, UDim>& Pu) {
  // B = B*Pu
  linearApproximation.dfdu =
      linearApproximation.dfdu * Pu;  // 临时矩阵不可避免。
}
