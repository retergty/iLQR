/******************************************************************************
Copyright (c) 2021, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

 * Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

/**
 * @file ChangeOfInputVariables.hpp
 * @brief 输入变量变换：对二次近似做 δu = Pu * δũ，将输入维度由 m 变为 p（Pu 为 m×p 矩阵）。
 */
#pragma once

#include "Types.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"

/**
 * @brief 对二次近似施加输入变换 δu = Pu * δũ，变换后状态维 n、输入维 p；Pu 为空矩阵时表示零输入。
 * @param [in,out] quadraticApproximation 待变换的二次近似（原地修改）。
 * @param [in] Pu 定义 δũ 范围的矩阵（m×p）。
 */
template <typename Scalar, int XDim, int UDim>
void changeOfInputVariables(ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> &quadraticApproximation,
                            const Matrix<Scalar, UDim, UDim> &Pu)
{
  // P = Pu'*P
  quadraticApproximation.dfdux = Pu.transpose() * quadraticApproximation.dfdux;

  // R = Pu' * R * Pu
  quadraticApproximation.dfduu = Pu.transpose() * quadraticApproximation.dfduu * Pu;

  // r = Pu' * r
  quadraticApproximation.dfdu = Pu.transpose() * quadraticApproximation.dfdu;
}

/** Applies the change of input variables to a linear system */
template <typename Scalar, int FDimisions, int XDim, int UDim>
void changeOfInputVariables(VectorFunctionLinearApproximation<Scalar, FDimisions, XDim, UDim> &linearApproximation, const Matrix<Scalar, UDim, UDim> &Pu)
{
  // B = B*Pu
  linearApproximation.dfdu = linearApproximation.dfdu * Pu; // temporary matrix unavoidable
}
