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
#include "StateAugmentedLagrangian.hpp"
#include <array>

/**
 * State Augmented Lagrangian penalty class combining a collection of constraint terms.
 *
 * This class collects a variable number of Augmented Lagrangian penalty terms and provides methods to get the
 * summed values and quadratic approximations. Each term can be accessed through its string name.
 */
template<typename Scalar, int XDimisions, int StateAugmentLagrangianNumbers>
class StateAugmentedLagrangianCollection
{
public:
  StateAugmentedLagrangianCollection() = default;

  /** Get state constraints and their penalties for each active term */
  std::array<LagrangianMetrics<Scalar>, StateAugmentLagrangianNumbers> getValue(const Scalar time, const Vector<Scalar, XDimisions>& state, const std::array<Multiplier<Scalar>, StateAugmentLagrangianNumbers>& termsMultiplier) const
  {
    std::array<LagrangianMetrics<Scalar>, StateAugmentLagrangianNumbers> termsConstraintPenalty;

    // int i = 0;
    // for (auto it = list_.begin();it != list_.end();++it)
    // {
    //   assert(i < StateAugmentLagrangianNumbers);

    //   termsConstraintPenalty[i] = it->getValue(time, state, termsMultiplier[i]);
    //   i++;
    // }

    for (int i = 0;i < num_;++i)
    {
      termsConstraintPenalty[i] = lagrangian_[i].getValue(time, state, termsMultiplier[i]);
    }
    return termsConstraintPenalty;
  }
  /** Get the sum of state Lagrangian penalties quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDimisions>& state, const std::array<Multiplier<Scalar>, StateAugmentLagrangianNumbers>& termsMultiplier) const
  {
    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> penalty;
    penalty.setZero();

    // // accumulate terms
    // int i = 0;
    // for (auto it = list_.begin();it != list_.end();++it)
    // {
    //   assert(i < StateAugmentLagrangianNumbers);

    //   penalty += it->getQuadraticApproximation(time, state, termsMultiplier[i]);
    //   i++;
    // }

    for (int i = 0;i < num_;++i)
    {
      penalty += lagrangian_[i].getQuadraticApproximation(time, state, termsMultiplier[i]);
    }
    return penalty;
  }

  /** Update Lagrange/penalty multipliers, and the penalty value for each active term. */
  void updateLagrangian(Scalar time, const Vector<Scalar, XDimisions>& state, std::array<LagrangianMetrics<Scalar>, StateAugmentLagrangianNumbers>& termsMetrics, std::array<Multiplier<Scalar>, StateAugmentLagrangianNumbers>& termsMultiplier) const
  {
    // int i = 0;
    // for (auto it = list_.begin();it != list_.end();++it)
    // {
    //   assert(i < StateAugmentLagrangianNumbers);

    //   // Multiplier<Scalar> updatedLagrangian;
    //   // std::tie(updatedLagrangian, termsMetrics[i].penalty) = it->updateLagrangian(time, state, termsMetrics[i].constraint, termsMultiplier[i]);
    //   // termsMultiplier[i] = updatedLagrangian;

    //   std::tie(termsMultiplier[i], termsMetrics[i].penalty) = it->updateLagrangian(time, state, termsMetrics[i].constraint, termsMultiplier[i]);

    //   i++;
    // }

    for (int i = 0;i < num_;++i)
    {
      std::tie(termsMultiplier[i], termsMetrics[i].penalty) = lagrangian_[i].updateLagrangian(time, state, termsMetrics[i].constraint, termsMultiplier[i]);
    }
  }

  /** Initialize Lagrange/penalty multipliers for each active term. */
  void initializeLagrangian(Scalar time, std::array<Multiplier<Scalar>, StateAugmentLagrangianNumbers>& termsMultiplier) const
  {
    //int i = 0;
    // for (auto it = list_.begin();it != list_.end();++it)
    // {
    //   assert(i < StateAugmentLagrangianNumbers);

    //   termsMultiplier[i] = it->initializeLagrangian(time);
    //   i++;
    // }
    for (int i = 0;i < num_; ++i)
    {
      termsMultiplier[i] = lagrangian_[i].initializeLagrangian(time);
    }
  }

  // add cost to list end
  void add(const StateAugmentedLagrangian<Scalar, XDimisions>& state_augment_lagrangian)
  {
    //list_.insert(list_.end(), state_augment_lagrangian);
    assert(num_ < StateAugmentLagrangianNumbers);
    lagrangian_[num_] = state_augment_lagrangian;
    num_++;
  }

private:
  int num_{ 0 };
  //IntrusiveList<StateAugmentedLagrangian<Scalar, XDimisions>> list_;
  std::array<StateAugmentedLagrangian<Scalar, XDimisions>, StateAugmentLagrangianNumbers> lagrangian_;
};
