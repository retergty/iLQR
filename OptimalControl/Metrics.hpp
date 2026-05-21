/**
 * @file Metrics.hpp
 * @brief 单时刻的指标：代价、动力学违反、各约束的拉格朗日项（惩罚与约束值）。
 */
#pragma once

#include <cmath>

#include "LagrangianMetrics.hpp"
#include "Types.hpp"

/**
 * @brief
 * 单时刻的代价、动力学违反与各类约束的拉格朗日指标（等式/不等式、状态/状态-输入）。
 * @tparam Scalar 标量类型。
 * @tparam Dims 状态/输入维度配置。
 * @tparam Layout 单时刻约束布局。
 */
template <typename Scalar, typename Dims, typename Layout>
struct Metrics {
  static constexpr int XDim = Dims::XDim;
  static constexpr int StateEq = Layout::StateEq;
  static constexpr int StateIneq = Layout::StateIneq;
  static constexpr int StateInputEq = Layout::StateInputEq;
  static constexpr int StateInputIneq = Layout::StateInputIneq;

  /** @brief 该时刻总代价。 */
  Scalar cost;

  /** @brief 动力学违反向量。 */
  Vector<Scalar, XDim> dynamicsViolation;

  /** @brief 状态等式约束的拉格朗日项数组。 */
  std::array<LagrangianMetrics<Scalar>, StateEq> stateEqLagrangian;
  /** @brief 状态不等式约束的拉格朗日项数组。 */
  std::array<LagrangianMetrics<Scalar>, StateIneq> stateIneqLagrangian;
  /** @brief 状态-输入等式约束的拉格朗日项数组。 */
  std::array<LagrangianMetrics<Scalar>, StateInputEq>
      stateInputEqLagrangian;
  /** @brief 状态-输入不等式约束的拉格朗日项数组。 */
  std::array<LagrangianMetrics<Scalar>, StateInputIneq>
      stateInputIneqLagrangian;

  /** @brief 与另一 Metrics 交换各成员。 */
  void swap(Metrics& other) {
    std::swap(cost, other.cost);
    dynamicsViolation.swap(other.dynamicsViolation);
    stateEqLagrangian.swap(other.stateEqLagrangian);
    stateIneqLagrangian.swap(other.stateIneqLagrangian);
    stateInputEqLagrangian.swap(other.stateInputEqLagrangian);
    stateInputIneqLagrangian.swap(other.stateInputIneqLagrangian);
  }

  /** @brief 清空/重置：代价置零、动力学违反置零；拉格朗日数组由
   * computeRolloutMetrics 覆盖。 */
  void clear() {
    cost = 0;
    dynamicsViolation.setZero();
  }
};

/**
 * @brief 求一组 LagrangianMetrics 的惩罚值之和。
 * @param [in] metricsArray 拉格朗日指标数组。
 * @return 所有元素的 penalty 之和。
 */
template <typename Scalar, size_t ArrayLen>
inline Scalar sumPenalties(
    const std::array<LagrangianMetrics<Scalar>, ArrayLen>& metricsArray) {
  Scalar s = 0.0;
  // std::for_each(metricsArray.begin(), metricsArray.end(), [&s](const
  // LagrangianMetrics<Scalar>& m) { s += m.penalty; });
  for (size_t i = 0; i < ArrayLen; ++i) {
    s += metricsArray[i].penalty;
  }
  return s;
}

/**
 * @brief 求一组 LagrangianMetrics 的约束值绝对值之和（标量约束）。
 * @param [in] metricsArray 拉格朗日指标数组。
 * @return 各元素 constraint 的绝对值之和。
 */
template <typename Scalar, size_t ArrayLen>
inline Scalar constraintsSquaredNorm(
    const std::array<LagrangianMetrics<Scalar>, ArrayLen>& metricsArray) {
  Scalar s = 0.0;
  // std::for_each(metricsArray.begin(), metricsArray.end(), [&s](const
  // LagrangianMetrics& m) { s += m.constraint.squaredNorm(); });
  for (int i = 0; i < ArrayLen; ++i) {
    s += std::abs(metricsArray[i].constraint);
  }
  return s;
}

/**
 * @brief 求等式约束违反向量的平方范数（SSE）。
 * @param [in] eqConstraint 等式约束违反向量。
 * @return eqConstraint.squaredNorm()，Dimisions==0 时为 0。
 */
template <typename Scalar, int Dimisions>
inline Scalar getEqConstraintsSSE(
    const Vector<Scalar, Dimisions>& eqConstraint) {
  static_assert(Dimisions >= 0);

  if constexpr (Dimisions == 0) {
    return 0;
  } else {
    return eqConstraint.squaredNorm();
  }
}