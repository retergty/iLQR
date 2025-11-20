/******************************************************************************
Copyright (c) 2020, Farbod Farshidian. All rights reserved.

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
 * The struct contains Riccati equation modification terms.
 */
template <typename Scalar, int XDimisions, int UDimisions>
struct RiccatiModification
{
  Scalar time_ = 0.0;

  Matrix<Scalar, XDimisions, XDimisions> deltaQm_;
  // Matrix<Scalar, XDimisions, UDimisions> deltaGm_;
  // Vector<Scalar, UDimisions> deltaGv_;

  /** The Hessian matrix of the Hamiltonian, \f$Hm\f$. */
  Matrix<Scalar, UDimisions, UDimisions> hamiltonianHessian_;
  // /** The right pseudo-inverse of \f$Dm\f$ */
  // Matrix<Scalar, UDimisions, XDimisions> constraintRangeProjector_;
  // /** \f$DmNull inv(DmNull^T * Hm * DmNull) * DmNull^T = (I - invHm * Dm^T * inv(Dm * invHm * Dm^T) * Dm) * invHm\f$ */
  Matrix<Scalar, UDimisions, UDimisions> constraintNullProjector_;
};
