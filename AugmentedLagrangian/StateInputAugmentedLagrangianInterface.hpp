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
#include "Multiplier.hpp"
#include "LagrangianMetrics.hpp"
#include "QuadraticApproximation.hpp"
#include "IntrusiveList.hpp"

/** The base class for Augmented Lagrangian penalty of state-input constraint. */
template<typename Scalar, int XDimisions, int UDimisions>
class StateInputAugmentedLagrangianInterface : IntrusiveListNode<StateInputAugmentedLagrangianInterface<Scalar, XDimisions, UDimisions>>
{
public:
  StateInputAugmentedLagrangianInterface() = default;
  virtual ~StateInputAugmentedLagrangianInterface() = default;

  /** Get the constraint and its penalty value */
  virtual LagrangianMetrics<Scalar> getValue(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input, const Multiplier<Scalar>& lagrangian) const = 0;

  /** Get the constraint's penalty quadratic approximation */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const Multiplier<Scalar>& lagrangian) const = 0;

  /** Update Lagrange/penalty multipliers and the penalty function value. */
  virtual std::pair<Multiplier<Scalar>, Scalar> updateLagrangian(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const Scalar constraint, const Multiplier<Scalar>& lagrangian) const = 0;

  /** Initialize Lagrange/penalty multipliers. */
  virtual Multiplier<Scalar> initializeLagrangian(Scalar time) const = 0;
};
