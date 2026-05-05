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
 * @file StateInputAugmentedLagrangianCollection.hpp
 * @brief 状态-输入增广拉格朗日集合：汇总多个状态-输入约束的惩罚项，提供总取值与总二次近似。
 */
#pragma once

#include "Types.hpp"
#include "StateInputAugmentedLagrangian.hpp"
#include <array>

/**
 * @brief 状态-输入增广拉格朗日惩罚项集合：对多个 StateInputAugmentedLagrangian 求和。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @tparam StateInputAugmentLagrangianNumbers 项数。
 */
template<typename Scalar, int XDim, int UDim, int StateInputAugmentLagrangianNumbers>
class StateInputAugmentedLagrangianCollection
{
public:
  StateInputAugmentedLagrangianCollection() = default;

  /** @brief 获取各激活项的状态-输入约束与惩罚值数组。 */
  std::array<LagrangianMetrics<Scalar>, StateInputAugmentLagrangianNumbers> getValue(const Scalar time, const Vector<Scalar, XDim>& state, const Vector<Scalar, UDim>& input, const std::array<Multiplier<Scalar>, StateInputAugmentLagrangianNumbers>& termsMultiplier) const
  {
    std::array<LagrangianMetrics<Scalar>, StateInputAugmentLagrangianNumbers> termsConstraintPenalty;

    // int i = 0;
    // for (auto it = list_.begin();it != list_.end();++it)
    // {
    //   assert(i < StateInputAugmentLagrangianNumbers);

    //   termsConstraintPenalty[i] = it->getValue(time, state, input, termsMultiplier[i]);
    //   i++;
    // }
    // return termsConstraintPenalty;

    for (int i = 0;i < num_;++i)
    {
      termsConstraintPenalty[i] = lagrangian_[i].getValue(time, state, termsMultiplier[i]);
    }
    return termsConstraintPenalty;
  }

  /** Get the sum of state-input Lagrangian penalties quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDim>& state, const Vector<Scalar, UDim>& input,
    const std::array<Multiplier<Scalar>, StateInputAugmentLagrangianNumbers>& termsMultiplier) const
  {
    ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> penalty;
    penalty.setZero();

    // // accumulate terms
    // int i = 0;
    // for (auto it = list_.begin();it != list_.end();++it)
    // {
    //   assert(i < StateInputAugmentLagrangianNumbers);

    //   penalty += it->getQuadraticApproximation(time, state, input, termsMultiplier[i]);
    //   i++;
    // }

    for (int i = 0;i < num_;++i)
    {
      penalty += lagrangian_[i].getQuadraticApproximation(time, state, termsMultiplier[i]);
    }
    return penalty;
  }

  /** Update Lagrange/penalty multipliers and the penalty value for each active term. */
  void updateLagrangian(const Scalar time, const Vector<Scalar, XDim>& state, const Vector<Scalar, UDim>& input, std::array<LagrangianMetrics<Scalar>, StateInputAugmentLagrangianNumbers>& termsMetrics,
    std::array<Multiplier<Scalar>, StateInputAugmentLagrangianNumbers>& termsMultiplier) const
  {
    // int i = 0;
    // for (auto it = list_.begin();it != list_.end();++it)
    // {
    //   assert(i < StateInputAugmentLagrangianNumbers);

    //   // Multiplier<Scalar> updatedLagrangian;
    //   // std::tie(updatedLagrangian, termsMetrics[i].penalty) = it->updateLagrangian(time, state, termsMetrics[i].constraint, termsMultiplier[i]);
    //   // termsMultiplier[i] = updatedLagrangian;

    //   std::tie(termsMultiplier[i], termsMetrics[i].penalty) = it->updateLagrangian(time, state, input, termsMetrics[i].constraint, termsMultiplier[i]);

    //   i++;
    // }
    for (int i = 0;i < num_;++i)
    {
      std::tie(termsMultiplier[i], termsMetrics[i].penalty) = lagrangian_[i].updateLagrangian(time, state, input, termsMetrics[i].constraint, termsMultiplier[i]);
    }
  }

  /** Initialize Lagrange/penalty multipliers for each active term. */
  void initializeLagrangian(const Scalar time, std::array<Multiplier<Scalar>, StateInputAugmentLagrangianNumbers>& termsMultiplier) const
  {
    // int i = 0;
    // for (auto it = list_.begin();it != list_.end();++it)
    // {
    //   assert(i < StateInputAugmentLagrangianNumbers);

    //   termsMultiplier[i] = it->initializeLagrangian(time);
    //   i++;
    // }

    for (int i = 0;i < num_; ++i)
    {
      termsMultiplier[i] = lagrangian_[i].initializeLagrangian(time);
    }
  }

  // add cost to list end
  void add(const StateInputAugmentedLagrangian<Scalar, XDim, UDim>& state_input_augment_lagrangian)
  {
    // list_.insert(list_.end(), state_input_augment_lagrangian);
    // num_++;

    assert(num_ < StateInputAugmentLagrangianNumbers);
    lagrangian_[num_] = state_input_augment_lagrangian;
    num_++;
  }

private:
  int num_{ 0 };
  //IntrusiveList<StateInputAugmentedLagrangian<Scalar, XDim, UDim>> list_;
  std::array<StateInputAugmentedLagrangian<Scalar, XDim, UDim>, StateInputAugmentLagrangianNumbers> lagrangian_;
};