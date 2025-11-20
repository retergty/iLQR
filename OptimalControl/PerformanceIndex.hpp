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
#include "Numerics.hpp"
#include "Metrics.hpp"

/**
 * Defines the performance indices for a rollout
 */
template<typename Scalar>
struct PerformanceIndex {
  /** The merit function of a rollout. */
  Scalar merit = 0.0;

  /** The total cost of a rollout. */
  Scalar cost = 0.0;

  /** Sum of Squared Error (SSE) of the dual feasibilities:
   * - Final: squared norm of violation in the dual feasibilities
   * - PreJumps: sum of squared norm of violation in the dual feasibilities
   * - Intermediates: sum of squared norm of violation in the dual feasibilities
   */
  Scalar dualFeasibilitiesSSE = 0.0;

  /** Sum of Squared Error (SSE) of system dynamics violation */
  Scalar dynamicsViolationSSE = 0.0;

  /** Sum of equality Lagrangians:
   * - Final: penalty for violation in state equality constraints
   * - PreJumps: penalty for violation in state equality constraints
   * - Intermediates: penalty for violation in state/state-input equality constraints
   */
  Scalar equalityLagrangian = 0.0;

  /** Sum of inequality Lagrangians:
   * - Final: penalty for violation in state inequality constraints
   * - PreJumps: penalty for violation in state inequality constraints
   * - Intermediates: penalty for violation in state/state-input inequality constraints
   */
  Scalar inequalityLagrangian = 0.0;

  /** Add performance indices. */
  PerformanceIndex& operator+=(const PerformanceIndex& rhs)
  {
    this->merit += rhs.merit;
    this->cost += rhs.cost;
    this->dualFeasibilitiesSSE += rhs.dualFeasibilitiesSSE;
    this->dynamicsViolationSSE += rhs.dynamicsViolationSSE;
    this->equalityLagrangian += rhs.equalityLagrangian;
    this->inequalityLagrangian += rhs.inequalityLagrangian;
    return *this;
  }

  /** Multiply by a scalar. */
  PerformanceIndex& operator*=(const Scalar c)
  {
    this->merit *= c;
    this->cost *= c;
    this->dualFeasibilitiesSSE *= c;
    this->dynamicsViolationSSE *= c;
    this->equalityLagrangian *= c;
    this->inequalityLagrangian *= c;
    return *this;
  }

  /** Returns true if *this is approximately equal to other, within the precision determined by prec. */
  bool isApprox(const PerformanceIndex& other, const Scalar prec = 1e-8) const
  {
    return numerics::almost_eq(this->merit, other.merit, prec) && numerics::almost_eq(this->cost, other.cost, prec) &&
      numerics::almost_eq(this->dualFeasibilitiesSSE, other.dualFeasibilitiesSSE, prec) &&
      numerics::almost_eq(this->dynamicsViolationSSE, other.dynamicsViolationSSE, prec) &&
      numerics::almost_eq(this->equalityLagrangian, other.equalityLagrangian, prec) &&
      numerics::almost_eq(this->inequalityLagrangian, other.inequalityLagrangian, prec);
  }
};

/** Add performance indices. */
template<typename Scalar>
inline PerformanceIndex<Scalar> operator+(PerformanceIndex<Scalar> lhs, const PerformanceIndex<Scalar>& rhs) {
  lhs += rhs;  // copied lhs, add rhs to it.
  return lhs;
}

/** Multiply by a scalar. */
template<typename Scalar>
inline PerformanceIndex<Scalar> operator*(PerformanceIndex<Scalar> lhs, const Scalar c) {
  lhs *= c;  // copied lhs
  return lhs;
}

/** Multiply by a scalar. */
template<typename Scalar>
inline PerformanceIndex<Scalar> operator*(const Scalar c, PerformanceIndex<Scalar> rhs) {
  rhs *= c;  // copied rhs
  return rhs;
}

/** Swaps performance indices. */
template<typename Scalar>
void swap(PerformanceIndex<Scalar>& lhs, PerformanceIndex<Scalar>& rhs)
{
  std::swap(lhs.merit, rhs.merit);
  std::swap(lhs.cost, rhs.cost);
  std::swap(lhs.dualFeasibilitiesSSE, rhs.dualFeasibilitiesSSE);
  std::swap(lhs.dynamicsViolationSSE, rhs.dynamicsViolationSSE);
  std::swap(lhs.equalityLagrangian, rhs.equalityLagrangian);
  std::swap(lhs.inequalityLagrangian, rhs.inequalityLagrangian);
}

/** Computes the PerformanceIndex based on a given continuous-time Metrics. */
template <typename Scalar, int XDimisions, int UDimisions, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains>
PerformanceIndex<Scalar> toPerformanceIndex(const Metrics<Scalar, XDimisions, UDimisions, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>& m)
{
  PerformanceIndex<Scalar> performanceIndex;
  performanceIndex.merit = 0.0;  // left for the solver to fill
  performanceIndex.cost = m.cost;
  performanceIndex.dualFeasibilitiesSSE = 0.0;  // left for the solver to fill
  performanceIndex.dynamicsViolationSSE = getEqConstraintsSSE(m.dynamicsViolation);
  performanceIndex.equalityLagrangian = sumPenalties(m.stateEqLagrangian) + sumPenalties(m.stateInputEqLagrangian);
  performanceIndex.inequalityLagrangian = sumPenalties(m.stateIneqLagrangian) + sumPenalties(m.stateInputIneqLagrangian);
  return performanceIndex;
}

/** Computes the PerformanceIndex based on a given discrete-time Metrics. */
template <typename Scalar, int XDimisions, int UDimisions, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains>
PerformanceIndex<Scalar> toPerformanceIndex(const Metrics<Scalar, XDimisions, UDimisions, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>& m, const Scalar dt)
{
  auto performanceIndex = toPerformanceIndex(m);
  //  performanceIndex.cost *= dt  no need since it is already considered in multiple_shooting::computeIntermediateMetrics()
  performanceIndex.dualFeasibilitiesSSE *= dt;
  performanceIndex.dynamicsViolationSSE *= dt;
  performanceIndex.inequalityConstraintsSSE *= dt;
  return performanceIndex;
}
