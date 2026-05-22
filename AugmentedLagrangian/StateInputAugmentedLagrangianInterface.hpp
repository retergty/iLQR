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
 * @file StateInputAugmentedLagrangianInterface.hpp
 * @brief 状态-输入增广拉格朗日接口：约束惩罚取值、二次近似与乘子初始化/更新。
 */
#pragma once

#include "LagrangianMetrics.hpp"
#include "Multiplier.hpp"
#include "QuadraticApproximation.hpp"
#include "Types.hpp"

/**
 * @brief
 * 状态-输入约束的增广拉格朗日惩罚接口：提供取值、二次近似、乘子更新与初始化。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int XDim, int UDim, int CDim>
class StateInputAugmentedLagrangianInterface {
 public:
  StateInputAugmentedLagrangianInterface() = default;
  virtual ~StateInputAugmentedLagrangianInterface() = default;

  /** Get the constraint and its penalty value */
  virtual LagrangianMetrics<Scalar, CDim> getValue(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const Multiplier<Scalar, CDim>& lagrangian) const = 0;

  /** Get the constraint's penalty quadratic approximation */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const Multiplier<Scalar, CDim>& lagrangian) const = 0;

  /** Update Lagrange/penalty multipliers and the penalty function value. */
  virtual std::pair<Multiplier<Scalar, CDim>, Scalar> updateLagrangian(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input, const Vector<Scalar, CDim>& constraint,
      const Multiplier<Scalar, CDim>& lagrangian) const = 0;

  /** Initialize Lagrange/penalty multipliers. */
  virtual Multiplier<Scalar, CDim> initializeLagrangian(Scalar time) const = 0;
};
