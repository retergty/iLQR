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
#include "StateConstraint.hpp"
#include "Types.hpp"

/**
 * Linear state-only constraint
 */
template <typename Scalar, int XDim, int CDim>
class LinearStateConstraint : public StateConstraint<Scalar, XDim, CDim> {
 public:
  /**
   * Constructor
   *
   * @param [in] h: Constant term in F * x + h = 0
   * @param [in] F: x factor in F * x + h = 0
   */
  LinearStateConstraint(const Vector<Scalar, CDim>& h,
                        const Matrix<Scalar, CDim, XDim>& F)
      : StateConstraint<Scalar, XDim, CDim>(ConstraintOrder::Linear),
        h_(h),
        F_(F) {}

  ~LinearStateConstraint() override = default;

  Vector<Scalar, CDim> getValue(const Scalar time,
                                const Vector<Scalar, XDim>& state) const final {
    (void)time;
    Vector<Scalar, CDim> g = h_;
    g += F_ * state;
    return g;
  }

  VectorFunctionLinearApproximation<Scalar, CDim, XDim, 0>
  getLinearApproximation(const Scalar time,
                         const Vector<Scalar, XDim>& state) const final {
    (void)time;
    VectorFunctionLinearApproximation<Scalar, CDim, XDim, 0> g;
    g.f = h_;
    g.f += F_ * state;
    g.dfdx = F_;
    return g;
  }

 public:
  Vector<Scalar, CDim> h_; /**< State only constraint */
  Matrix<Scalar, CDim, XDim>
      F_; /**< State only constraint derivative wrt. state */
};
