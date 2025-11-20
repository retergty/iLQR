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
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"

/**
 * The penalty function interface class is used to penalize constraint violation by adding a penalty term to the cost function.
 * We assume that the penalty function is convex. In general, the penalty is a function of time, Lagrange multiplier, and
 * constraint violation.
 */
template<typename Scalar>
class AugmentedPenaltyBase {
public:
  /** Default constructor */
  AugmentedPenaltyBase() = default;

  /** Default destructor */
  virtual ~AugmentedPenaltyBase() = default;

  /**
   * Compute the penalty value at a certain constraint value.
   *
   * @param [in] t: The time that the constraint is evaluated.
   * @param [in] l: The Lagrange multiplier.
   * @param [in] h: Constraint value.
   * @return penalty cost.
   */
  virtual Scalar getValue(const Scalar t, const Scalar l, const Scalar h) const = 0;
  
  /**
   * Compute the penalty derivative at a certain constraint value.
   *
   * @param [in] t: The time that the constraint is evaluated.
   * @param [in] l: The Lagrange multiplier.
   * @param [in] h: Constraint value.
   * @return penalty derivative with respect to constraint value.
   */
  virtual Scalar getDerivative(const Scalar t, const Scalar l, const Scalar h) const = 0;

  /**
   * Compute the penalty second derivative at a certain constraint value.
   *
   * @param [in] t: The time that the constraint is evaluated.
   * @param [in] l: The Lagrange multiplier.
   * @param [in] h: Constraint value.
   * @return penalty second derivative with respect to constraint value.
   */
  virtual Scalar getSecondDerivative(const Scalar t, const Scalar l, const Scalar h) const = 0;

  /**
   * Updates the Lagrange multiplier.
   *
   * @param [in] t: The time stamp.
   * @param [in] l: The Lagrange multiplier.
   * @param [in] h: Constraint values.
   * @return updated Lagrange multiplier.
   */
  virtual Scalar updateMultiplier(const Scalar t, const Scalar l, const Scalar h) const = 0;

  /**
   * Initializes the Lagrange multiplier.
   *
   * @return initial value of the Lagrange multiplier.
   */
  virtual Scalar initializeMultiplier() const = 0;

  AugmentedPenaltyBase(const AugmentedPenaltyBase& other) = default;
};
