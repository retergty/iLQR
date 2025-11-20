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

#pragma once

#include "Types.hpp"

/**
 * @brief The Hessian matrix correction strategy
 * Enum used in selecting either DIAGONAL_SHIFT, CHOLESKY_MODIFICATION, EIGENVALUE_MODIFICATION, or GERSHGORIN_MODIFICATION strategies.
 */
enum class HessianCorrectionStrategy
{
  DIAGONAL_SHIFT
  // CHOLESKY_MODIFICATION,
  // EIGENVALUE_MODIFICATION,
  // GERSHGORIN_MODIFICATION
};

/**
 * Shifts the Hessian based on the strategy defined by Line_Search::hessianCorrectionStrategy_.
 *
 * @param [in] strategy: Hessian matrix correction strategy.
 * @param [in, out] matrix: The Hessian matrix.
 * @param [in] minEigenvalue: The minimum expected eigenvalue after correction.
 */
template <typename Scalar, int Dimisions>
void shiftHessian(HessianCorrectionStrategy strategy, Matrix<Scalar, Dimisions, Dimisions> &matrix, Scalar minEigenvalue = 1e-6)
{
  switch (strategy)
  {
  case HessianCorrectionStrategy::DIAGONAL_SHIFT:
  {
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
