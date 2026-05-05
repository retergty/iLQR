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

/**
 * @file PerformanceIndex.hpp
 * @brief 单次 rollout 的性能指标：代价、merit、对偶可行性 SSE、动力学违反、等式/不等式拉格朗日等。
 */
#pragma once

#include "Types.hpp"
#include "Numerics.hpp"
#include "Metrics.hpp"

/**
 * @brief 单次 rollout 的性能指标汇总，用于收敛判断与 merit 计算。
 * @tparam Scalar 标量类型。
 */
template<typename Scalar>
struct PerformanceIndex {
  /** @brief 该次 rollout 的 merit 函数值（用于线搜索与收敛判据）。 */
  Scalar merit = 0.0;

  /** @brief 该次 rollout 的总代价。 */
  Scalar cost = 0.0;

  /** @brief 对偶可行性的误差平方和（终端/中间时刻违反的平方范数之和）。 */
  Scalar dualFeasibilitiesSSE = 0.0;

  /** @brief 系统动力学违反的误差平方和。 */
  Scalar dynamicsViolationSSE = 0.0;

  /** @brief 等式约束拉格朗日项之和（状态等式、状态-输入等式惩罚）。 */
  Scalar equalityLagrangian = 0.0;

  /** @brief 不等式约束拉格朗日项之和（状态不等式、状态-输入不等式惩罚）。 */
  Scalar inequalityLagrangian = 0.0;

  /** @brief 将另一性能指标逐项加到本对象。 */
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

  /** @brief 将本性能指标各分量乘以标量 c。 */
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

  /**
   * @brief 判断与另一性能指标是否在给定精度内近似相等。
   * @param [in] other 另一性能指标。
   * @param [in] prec 数值比较精度，默认 1e-8。
   * @return 各分量均在 prec 内近似相等则返回 true。
   */
  bool isApprox(const PerformanceIndex& other, const Scalar prec = 1e-8) const
  {
    return numerics::almost_eq(this->merit, other.merit, prec) && numerics::almost_eq(this->cost, other.cost, prec) &&
      numerics::almost_eq(this->dualFeasibilitiesSSE, other.dualFeasibilitiesSSE, prec) &&
      numerics::almost_eq(this->dynamicsViolationSSE, other.dynamicsViolationSSE, prec) &&
      numerics::almost_eq(this->equalityLagrangian, other.equalityLagrangian, prec) &&
      numerics::almost_eq(this->inequalityLagrangian, other.inequalityLagrangian, prec);
  }
};

/** @brief 两个性能指标逐项相加，返回新对象。 */
template<typename Scalar>
inline PerformanceIndex<Scalar> operator+(PerformanceIndex<Scalar> lhs, const PerformanceIndex<Scalar>& rhs) {
  lhs += rhs;
  return lhs;
}

/** @brief 性能指标各分量乘以标量 c。 */
template<typename Scalar>
inline PerformanceIndex<Scalar> operator*(PerformanceIndex<Scalar> lhs, const Scalar c) {
  lhs *= c;
  return lhs;
}

/** @brief 标量 c 乘以性能指标（与 operator*(lhs, c) 等价）。 */
template<typename Scalar>
inline PerformanceIndex<Scalar> operator*(const Scalar c, PerformanceIndex<Scalar> rhs) {
  rhs *= c;
  return rhs;
}

/** @brief 交换两个性能指标对象的内容。 */
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

/**
 * @brief 由单点连续时间 Metrics 计算对应的 PerformanceIndex（代价、动力学违反、等式/不等式拉格朗日等）。
 * @param [in] m 该时刻的 Metrics。
 * @return 填充 cost、dynamicsViolationSSE、equalityLagrangian、inequalityLagrangian 等的 PerformanceIndex。
 */
template <typename Scalar, int XDim, int UDim, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains>
PerformanceIndex<Scalar> toPerformanceIndex(const Metrics<Scalar, XDim, UDim, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>& m)
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

/**
 * @brief 由单点离散时间 Metrics 及时间步长 dt 计算 PerformanceIndex，部分 SSE 项乘以 dt。
 * @param [in] m 该时刻的 Metrics。
 * @param [in] dt 时间步长。
 * @return 乘以 dt 后的 PerformanceIndex（cost 已在 computeIntermediateMetrics 中考虑）。
 */
template <typename Scalar, int XDim, int UDim, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains>
PerformanceIndex<Scalar> toPerformanceIndex(const Metrics<Scalar, XDim, UDim, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>& m, const Scalar dt)
{
  auto performanceIndex = toPerformanceIndex(m);
  performanceIndex.dualFeasibilitiesSSE *= dt;
  performanceIndex.dynamicsViolationSSE *= dt;
  return performanceIndex;
}
