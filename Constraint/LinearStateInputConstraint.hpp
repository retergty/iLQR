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

#pragma once

#include "StateInputConstraint.hpp"
#include "Types.hpp"

/**
 * Linear state-input constraint
 */
template <typename Scalar, int XDim, int UDim, int CDim>
class LinearStateInputConstraint
    : public StateInputConstraint<Scalar, XDim, UDim, CDim> {
 public:
  /**
   * Constructor
   *
   * @param[in] e: Constant term in C * x + D * u + e = 0
   * @param[in] C: x factor in C * x + D * u + e = 0
   * @param[in] D: u factor in C * x + D * u + e = 0
   */
  LinearStateInputConstraint(const Vector<Scalar, CDim>& e,
                             const Matrix<Scalar, CDim, XDim>& C,
                             const Matrix<Scalar, CDim, UDim>& D)
      : StateInputConstraint<Scalar, XDim, UDim, CDim>(ConstraintOrder::Linear),
        e_(e),
        C_(C),
        D_(D) {}

  ~LinearStateInputConstraint() override = default;

  Vector<Scalar, CDim> getValue(const Scalar time,
                                const Vector<Scalar, XDim>& state,
                                const Vector<Scalar, UDim>& input) const final {
    (void)time;
    Vector<Scalar, CDim> g = e_;
    g += C_ * state;
    g += D_ * input;
    return g;
  }

  VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim>
  getLinearApproximation(const Scalar time, const Vector<Scalar, XDim>& state,
                         const Vector<Scalar, UDim>& input) const final {
    (void)time;
    VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim> g;
    g.f = e_;
    g.f += C_ * state;
    g.f += D_ * input;
    g.dfdx = C_;
    g.dfdu = D_;
    return g;
  }

 public:
  Vector<Scalar, CDim> e_; /**< State input constraint */
  Matrix<Scalar, CDim, XDim>
      C_; /**< State input constraint derivative wrt. state */
  Matrix<Scalar, CDim, UDim>
      D_; /**< State input constraint derivative wrt. input */
};
