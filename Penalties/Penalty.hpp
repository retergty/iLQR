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
#include "AugmentedPenaltyBase.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include <tuple>
/**
 *   A helper class that implements the penalty for multidimensional constraint
 *   \f$ h_i(x, u) \quad \forall  i \in [1,..,M] \f$
 *
 *   penalty(t, x, u) = \f$ \sum_{i=1}^{M} p(t, h_i(x, u)) \f$
 *
 *   This class uses the chain rule to compute the second-order approximation of the constraint-penalty. In the case that the
 *   second-order approximation of constraint is not provided, it employs a Gauss-Newton approximation technique which only
 *   relies on the first-order approximation. In general, the penalty function can be a function of time.
*/
template<typename Scalar, int XDimisions, int UDimisions>
class Penalty final
{
public:
  /**
   * Constructor
   */
  Penalty(AugmentedPenaltyBase<Scalar>* penaltyPtr) : penalty_ptr_(penaltyPtr) {};

  /** Default destructor */
  ~Penalty() = default;

  /** Copy constructor */
  Penalty(const Penalty& other) = delete;

  /**
   * Get the penalty cost.
   *
   * @param [in] t: The time that the constraint is evaluated.
   * @param [in] h: Vector of inequality constraint values.
   * @return Penalty: The penalty cost.
   */
  Scalar getValue(const Scalar t, const Scalar h, const Scalar l) const
  {
    return penalty_ptr_->getValue(t, l, h);
  }

  /**
   * Get the derivative of the penalty cost.
   * Implements the chain rule between the inequality constraint and penalty function.
   *
   * @param [in] t: The time that the constraint is evaluated.
   * @param [in] h: The constraint linear approximation.
   * @return The penalty cost quadratic approximation.
   */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> getQuadraticApproximation(const Scalar t, const ScalarFunctionLinearApproximation<Scalar, XDimisions, UDimisions>& h, const Scalar l) const
  {
    Scalar penaltyValue = 0.0;
    Scalar penaltyDerivative, penaltySecondDerivative;
    std::tie(penaltyValue, penaltyDerivative, penaltySecondDerivative) = getPenaltyValue1stDev2ndDev(t, h.f, l);
    const Vector<Scalar, XDimisions> penaltySecondDev_dhdx = penaltySecondDerivative * h.dfdx;

    // to make sure that dfdux in the state-only case has a right size
    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> penaltyApproximation;

    penaltyApproximation.f = penaltyValue;
    penaltyApproximation.dfdx = h.dfdx * penaltyDerivative;
    penaltyApproximation.dfdxx = h.dfdx * penaltySecondDev_dhdx.transpose();

    if constexpr (UDimisions > 0)
    {
      penaltyApproximation.dfdu = h.dfdu * penaltyDerivative;
      penaltyApproximation.dfdux = h.dfdu * penaltySecondDev_dhdx.transpose();
      penaltyApproximation.dfduu = h.dfdu * penaltySecondDerivative * h.dfdu.transpose();
    }

    return penaltyApproximation;
  }

  /**
   * Get the derivative of the penalty cost.
   * Implements the chain rule between the inequality constraint and penalty function.
   *
   * @param [in] t: The time that the constraint is evaluated.
   * @param [in] h: The constraint quadratic approximation.
   * @return The penalty cost quadratic approximation.
   */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> getQuadraticApproximation(Scalar t, const ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>& h,
    const Scalar l) const
  {
    const auto stateDim = h.dfdx.cols();
    const auto inputDim = h.dfdu.cols();
    const auto numConstraints = h.f.rows();

    Scalar penaltyValue = 0.0;
    Scalar penaltyDerivative, penaltySecondDerivative;
    std::tie(penaltyValue, penaltyDerivative, penaltySecondDerivative) = getPenaltyValue1stDev2ndDev(t, h.f, l);
    const Vector<Scalar, XDimisions> penaltySecondDev_dhdx = penaltySecondDerivative * h.dfdx;

    // to make sure that dfdux in the state-only case has a right size
    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> penaltyApproximation;

    penaltyApproximation.f = penaltyValue;
    penaltyApproximation.dfdx = h.dfdx * penaltyDerivative;
    penaltyApproximation.dfdxx = h.dfdx * penaltySecondDev_dhdx.transpose();

    penaltyApproximation.dfdxx += penaltyDerivative * h.dfdxx;

    if constexpr (UDimisions > 0) {
      penaltyApproximation.dfdu = h.dfdu * penaltyDerivative;
      penaltyApproximation.dfdux = h.dfdu * penaltySecondDev_dhdx.transpose();
      penaltyApproximation.dfduu = h.dfdu * penaltySecondDerivative * h.dfdu.transpose();

      penaltyApproximation.dfduu += penaltyDerivative * h.dfduu;
      penaltyApproximation.dfdux += penaltyDerivative * h.dfdux;
    }

    return penaltyApproximation;

  }

  /**
   * Updates the Lagrange multipliers.
   *
   * @param [in] t: The time stamp.
   * @param [in] l: The Lagrange multiplier.
   * @param [in] h: The constraint value.
   * @return updated Lagrange multipliers.
   */
  Scalar updateMultipliers(Scalar t, const Scalar h, const Scalar l) const
  {
    return penalty_ptr_->updateMultiplier(t, l, h);
  }

  /**
   * Initializes the Lagrange multipliers.
   * @return Initial Lagrange multipliers.
   */
  Scalar initializeMultipliers() const
  {
    return penalty_ptr_->initializeMultiplier();
  }

private:
  std::tuple<Scalar, Scalar, Scalar> getPenaltyValue1stDev2ndDev(const Scalar t, const Scalar h, const Scalar l) const
  {
    Scalar penaltyValue = penalty_ptr_->getValue(t, l, h);
    Scalar penaltyDerivative = penalty_ptr_->getDerivative(t, l, h);
    Scalar penaltySecondDerivative = penalty_ptr_->getSecondDerivative(t, l, h);

    return { penaltyValue, penaltyDerivative, penaltySecondDerivative };
  }

  AugmentedPenaltyBase<Scalar>* penalty_ptr_;
};
