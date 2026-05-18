/******************************************************************************
Copyright (c) 2017, Farbod Farshidian. All rights reserved.

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
 * @file HessianCorrection.hpp
 * @brief Hessian 修正策略：对角平移等，保证 Riccati 递推数值稳定性。
 */
#pragma once

#include "Types.hpp"

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
      matrix.diagonal().array() += minEigenvalue;
      break;
    }
  }
  // case HessianCorrectionStrategy::CHOLESKY_MODIFICATION: {
  //   LinearAlgebra::makePsdCholesky(matrix, minEigenvalue);
  //   break;
  // }
  // case HessianCorrectionStrategy::EIGENVALUE_MODIFICATION: {
  //   LinearAlgebra::makePsdEigenvalue(matrix, minEigenvalue);
  //   break;
  // }
  // case HessianCorrectionStrategy::GERSHGORIN_MODIFICATION: {
  //   LinearAlgebra::makePsdGershgorin(matrix, minEigenvalue);
  //   break;
  // }
}
