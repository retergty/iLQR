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
#include "QuadraticApproximation.hpp"
#include "LinearApproximation.hpp"

namespace qp_solver {

  /** Defines the quadratic cost and  linear dynamics at a give stage */
  template <typename Scalar, int XDim, int UDim>
  struct LinearQuadraticStage {
    using CostApproximation_t = ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>;
    using DynamicsApproximation_t = VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>;
    using ConstraintApproximation_t = VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, XDim, UDim>;
    /** Quadratic approximation of the cost */
    CostApproximation_t cost;
    /** Linear approximation of the dynamics */
    DynamicsApproximation_t dynamics;
    /** Linear approximation of the constraints */
    ConstraintApproximation_t constraints;

    LinearQuadraticStage()
    {
      constraints.f.resize(0);
      constraints.dfdx.resize(0, XDim);
      constraints.dfdu.resize(0, UDim);
    }

    LinearQuadraticStage(CostApproximation_t c, DynamicsApproximation_t d, ConstraintApproximation_t g)
      : cost(std::move(c)), dynamics(std::move(d)), constraints(std::move(g)) {
    }
  };

}  // namespace qp_solver

