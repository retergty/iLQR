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
 * @file Penalty.hpp
 * @brief 约束惩罚封装：基于 AugmentedPenaltyBase 计算惩罚值、二次近似及乘子更新。
 *
 * 对约束 h(x,u)，惩罚为 p(t, h, l)；本类用链式法则计算约束-惩罚的二阶近似，
 * 并委托底层惩罚类更新拉格朗日乘子。
 */
/**
 * @brief 单约束惩罚封装：取值、二次近似与乘子初始化/更新，均委托 penalty_ptr_。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template<typename Scalar, int XDim, int UDim>
class Penalty final
{
public:
  /** @brief 用增广惩罚基类指针构造。 */
  Penalty(AugmentedPenaltyBase<Scalar>* penaltyPtr) : penalty_ptr_(penaltyPtr) {};

  /** @brief 析构函数。 */
  ~Penalty() = default;

  /** @brief 禁止拷贝。 */
  Penalty(const Penalty& other) = delete;

  /**
   * @brief 获取惩罚代价值 p(t, h, l)。
   * @param [in] t 评估时间。
   * @param [in] h 约束值。
   * @param [in] l 拉格朗日乘子。
   * @return 惩罚值。
   */
  Scalar getValue(const Scalar t, const Scalar h, const Scalar l) const
  {
    return penalty_ptr_->getValue(t, l, h);
  }

  /**
   * @brief 由约束的线性近似经链式法则得到惩罚的二次近似。
   * @param [in] t 评估时间。
   * @param [in] h 约束的线性近似。
   * @param [in] l 拉格朗日乘子。
   * @return 惩罚的二次近似。
   */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> getQuadraticApproximation(const Scalar t, const ScalarFunctionLinearApproximation<Scalar, XDim, UDim>& h, const Scalar l) const
  {
    Scalar penaltyValue = 0.0;
    Scalar penaltyDerivative, penaltySecondDerivative;
    std::tie(penaltyValue, penaltyDerivative, penaltySecondDerivative) = getPenaltyValue1stDev2ndDev(t, h.f, l);
    const Vector<Scalar, XDim> penaltySecondDev_dhdx = penaltySecondDerivative * h.dfdx;

    // to make sure that dfdux in the state-only case has a right size
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> penaltyApproximation;

    penaltyApproximation.f = penaltyValue;
    penaltyApproximation.dfdx = h.dfdx * penaltyDerivative;
    penaltyApproximation.dfdxx = h.dfdx * penaltySecondDev_dhdx.transpose();

    if constexpr (UDim > 0)
    {
      penaltyApproximation.dfdu = h.dfdu * penaltyDerivative;
      penaltyApproximation.dfdux = h.dfdu * penaltySecondDev_dhdx.transpose();
      penaltyApproximation.dfduu = h.dfdu * penaltySecondDerivative * h.dfdu.transpose();
    }

    return penaltyApproximation;
  }

  /**
   * @brief 由约束的二次近似经链式法则得到惩罚的二次近似。
   * @param [in] t 评估时间。
   * @param [in] h 约束的二次近似。
   * @param [in] l 拉格朗日乘子。
   * @return 惩罚的二次近似。
   */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> getQuadraticApproximation(Scalar t, const ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& h,
    const Scalar l) const
  {
    const auto stateDim = h.dfdx.cols();
    const auto inputDim = h.dfdu.cols();
    const auto numConstraints = h.f.rows();

    Scalar penaltyValue = 0.0;
    Scalar penaltyDerivative, penaltySecondDerivative;
    std::tie(penaltyValue, penaltyDerivative, penaltySecondDerivative) = getPenaltyValue1stDev2ndDev(t, h.f, l);
    const Vector<Scalar, XDim> penaltySecondDev_dhdx = penaltySecondDerivative * h.dfdx;

    // to make sure that dfdux in the state-only case has a right size
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> penaltyApproximation;

    penaltyApproximation.f = penaltyValue;
    penaltyApproximation.dfdx = h.dfdx * penaltyDerivative;
    penaltyApproximation.dfdxx = h.dfdx * penaltySecondDev_dhdx.transpose();

    penaltyApproximation.dfdxx += penaltyDerivative * h.dfdxx;

    if constexpr (UDim > 0) {
      penaltyApproximation.dfdu = h.dfdu * penaltyDerivative;
      penaltyApproximation.dfdux = h.dfdu * penaltySecondDev_dhdx.transpose();
      penaltyApproximation.dfduu = h.dfdu * penaltySecondDerivative * h.dfdu.transpose();

      penaltyApproximation.dfduu += penaltyDerivative * h.dfduu;
      penaltyApproximation.dfdux += penaltyDerivative * h.dfdux;
    }

    return penaltyApproximation;

  }

  /**
   * @brief 根据约束值 h 与当前乘子 l 更新拉格朗日乘子。
   * @param [in] t 时间戳。
   * @param [in] h 约束值。
   * @param [in] l 当前拉格朗日乘子。
   * @return 更新后的乘子。
   */
  Scalar updateMultipliers(Scalar t, const Scalar h, const Scalar l) const
  {
    return penalty_ptr_->updateMultiplier(t, l, h);
  }

  /** @brief 初始化拉格朗日乘子。 @return 初始乘子。 */
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
