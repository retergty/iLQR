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

/**
 * @file StateAugmentedLagrangianCollection.hpp
 * @brief 仅状态增广拉格朗日集合：汇总多个仅状态约束的惩罚项，提供总取值与总二次近似。
 */
#pragma once

#include "Types.hpp"
#include "StateAugmentedLagrangian.hpp"
#include <array>

 /**
  * @brief 仅状态增广拉格朗日惩罚项集合：对多个 StateAugmentedLagrangian 求和。
  * @tparam Scalar 标量类型。
  * @tparam XDim 状态维度。
  * @tparam StateAugmentLagrangianNumbers 项数。
  */
template<typename Scalar, int XDim, int StateAugmentLagrangianNumbers>
class StateAugmentedLagrangianCollection
{
public:
  StateAugmentedLagrangianCollection() = default;

  /** @brief 获取各激活项的约束与惩罚值数组。 */
  std::array<LagrangianMetrics<Scalar>, StateAugmentLagrangianNumbers> getValue(const Scalar time, const Vector<Scalar, XDim>& state, const std::array<Multiplier<Scalar>, StateAugmentLagrangianNumbers>& termsMultiplier) const
  {
    std::array<LagrangianMetrics<Scalar>, StateAugmentLagrangianNumbers> termsConstraintPenalty;

    for (int i = 0;i < num_;++i)
    {
      termsConstraintPenalty[i] = lagrangian_[i]->getValue(time, state, termsMultiplier[i]);
    }
    return termsConstraintPenalty;
  }
  /** Get the sum of state Lagrangian penalties quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDim>& state, const std::array<Multiplier<Scalar>, StateAugmentLagrangianNumbers>& termsMultiplier) const
  {
    ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> penalty;
    penalty.setZero();

    for (int i = 0;i < num_;++i)
    {
      penalty += lagrangian_[i]->getQuadraticApproximation(time, state, termsMultiplier[i]);
    }
    return penalty;
  }

  /** Update Lagrange/penalty multipliers, and the penalty value for each active term. */
  void updateLagrangian(Scalar time, const Vector<Scalar, XDim>& state, std::array<LagrangianMetrics<Scalar>, StateAugmentLagrangianNumbers>& termsMetrics, std::array<Multiplier<Scalar>, StateAugmentLagrangianNumbers>& termsMultiplier) const
  {
    for (int i = 0;i < num_;++i)
    {
      std::tie(termsMultiplier[i], termsMetrics[i].penalty) = lagrangian_[i]->updateLagrangian(time, state, termsMetrics[i].constraint, termsMultiplier[i]);
    }
  }

  /** Initialize Lagrange/penalty multipliers for each active term. */
  void initializeLagrangian(Scalar time, std::array<Multiplier<Scalar>, StateAugmentLagrangianNumbers>& termsMultiplier) const
  {
    for (int i = 0;i < num_; ++i)
    {
      termsMultiplier[i] = lagrangian_[i]->initializeLagrangian(time);
    }
  }

  // add cost to list end
  void add(const StateAugmentedLagrangian<Scalar, XDim>* state_augment_lagrangian)
  {
    assert(num_ < StateAugmentLagrangianNumbers);
    lagrangian_[num_] = state_augment_lagrangian;
    num_++;
  }

private:
  int num_{ 0 };
  std::array<const StateAugmentedLagrangian<Scalar, XDim>*, StateAugmentLagrangianNumbers> lagrangian_;
};
