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
 * @file StateInputAugmentedLagrangian.hpp
 * @brief 状态-输入增广拉格朗日实现：绑定约束与惩罚，提供取值、二次近似与乘子更新。
 */
#pragma once

#include "Types.hpp"
#include "StateInputAugmentedLagrangianInterface.hpp"
#include "StateInputConstraint.hpp"
#include "AugmentedPenaltyBase.hpp"
#include "Penalty.hpp"
#include "LagrangianMetrics.hpp"
#include "QuadraticApproximation.hpp"

/**
 * @brief 状态-输入约束的增广拉格朗日惩罚实现：委托约束与 Penalty 计算取值与二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int XDim, int UDim>
class StateInputAugmentedLagrangian final : public StateInputAugmentedLagrangianInterface<Scalar, XDim, UDim>
{
public:
  /**
   * @brief 用约束指针与惩罚指针构造。
   * @param [in] constraintPtr 作为软约束的状态-输入约束。
   * @param [in] augmented_penalty 约束上的惩罚函数。
   */
  StateInputAugmentedLagrangian(StateInputConstraint<Scalar, XDim, UDim> *constraintPtr, AugmentedPenaltyBase<Scalar> *augmented_penalty) : constraint_ptr_(constraintPtr), penalty_(augmented_penalty) {};

  LagrangianMetrics<Scalar> getValue(const Scalar time, const Vector<Scalar, XDim> &state, const Vector<Scalar, UDim> &input, const Multiplier<Scalar> &multiplier) const override
  {
    const Scalar h = constraint_ptr_->getValue(time, state, input);
    const Scalar p = multiplier.penalty * penalty_.getValue(time, h, multiplier.lagrangian);
    return {p, h};
  }

  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDim> &state, const Vector<Scalar, UDim> &input, const Multiplier<Scalar> &multiplier) const override
  {
    switch (constraint_ptr_->getOrder())
    {
    case ConstraintOrder::Linear:
      return multiplier.penalty * penalty_.getQuadraticApproximation(time, constraint_ptr_->getLinearApproximation(time, state, input), multiplier.lagrangian);
    case ConstraintOrder::Quadratic:
      return multiplier.penalty * penalty_.getQuadraticApproximation(time, constraint_ptr_->getQuadraticApproximation(time, state, input), multiplier.lagrangian);
    default:
      return ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>();
    }
  }

  std::pair<Multiplier<Scalar>, Scalar> updateLagrangian(const Scalar time, const Vector<Scalar, XDim> &state, const Vector<Scalar, UDim> &input, const Scalar constraint, const Multiplier<Scalar> &multiplier) const override
  {
    (void)input;
    (void)state;

    const Multiplier<Scalar> updatedMultiplier{multiplier.penalty, penalty_.updateMultipliers(time, constraint, multiplier.lagrangian)};
    const Scalar penalty = updatedMultiplier.penalty * penalty_.getValue(time, constraint, updatedMultiplier.lagrangian);
    return {updatedMultiplier, penalty};
  }

  Multiplier<Scalar> initializeLagrangian(const Scalar time) const override
  {
    (void)time;
    return {1.0, penalty_.initializeMultipliers()};
  }

private:
  StateInputConstraint<Scalar, XDim, UDim>* constraint_ptr_;
  Penalty<Scalar, XDim, UDim> penalty_;
};
