/**
 * @file Metrics.hpp
 * @brief 单时刻的指标：代价、动力学违反、各约束的拉格朗日项（惩罚与约束值）。
 */
#pragma once

#include <tuple>

#include "OptimalControlData/LagrangianMetrics.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"
#include "iLQR/iLQRDescriptor.hpp"

template <typename Scalar, typename GroupLayout>
struct LagrangianMetricsGroup;

template <typename Scalar, typename... Terms>
struct LagrangianMetricsGroup<Scalar, ConstraintGroupLayout<Terms...>> {
  /** @brief 每个 constraint term 对应一个固定维度的拉格朗日指标。 */
  std::tuple<LagrangianMetrics<Scalar, Terms::CDim>...> terms;

  /** @brief 与另一 LagrangianMetricsGroup 交换内容。 */
  void swap(LagrangianMetricsGroup& other) { terms.swap(other.terms); }

  /** @brief 清空组内所有拉格朗日指标。 */
  void setZero() {
    std::apply([](auto&... metrics) { (metrics.setZero(), ...); }, terms);
  }
};

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
  using StateEqLayout = typename Layout::StateEqLayout;
  using StateIneqLayout = typename Layout::StateIneqLayout;
  using StateInputEqLayout = typename Layout::StateInputEqLayout;
  using StateInputIneqLayout = typename Layout::StateInputIneqLayout;

  static constexpr int StateEq = StateEqLayout::TotalDim;
  static constexpr int StateIneq = StateIneqLayout::TotalDim;
  static constexpr int StateInputEq = StateInputEqLayout::TotalDim;
  static constexpr int StateInputIneq = StateInputIneqLayout::TotalDim;

  /** @brief 默认构造为清空状态。 */
  Metrics() { clear(); }

  /** @brief 该时刻总代价。 */
  Scalar cost;

  /** @brief 动力学违反向量。 */
  Vector<Scalar, XDim> dynamicsViolation;

  /** @brief 状态等式约束的拉格朗日指标组。 */
  LagrangianMetricsGroup<Scalar, StateEqLayout> stateEqLagrangian;
  /** @brief 状态不等式约束的拉格朗日指标组。 */
  LagrangianMetricsGroup<Scalar, StateIneqLayout> stateIneqLagrangian;
  /** @brief 状态-输入等式约束的拉格朗日指标组。 */
  LagrangianMetricsGroup<Scalar, StateInputEqLayout> stateInputEqLagrangian;
  /** @brief 状态-输入不等式约束的拉格朗日指标组。 */
  LagrangianMetricsGroup<Scalar, StateInputIneqLayout> stateInputIneqLagrangian;

  /** @brief 与另一 Metrics 交换各成员。 */
  void swap(Metrics& other) {
    std::swap(cost, other.cost);
    dynamicsViolation.swap(other.dynamicsViolation);
    stateEqLagrangian.swap(other.stateEqLagrangian);
    stateIneqLagrangian.swap(other.stateIneqLagrangian);
    stateInputEqLagrangian.swap(other.stateInputEqLagrangian);
    stateInputIneqLagrangian.swap(other.stateInputIneqLagrangian);
  }

  /** @brief 清空/重置所有指标。 */
  void clear() {
    cost = 0;
    dynamicsViolation.setZero();
    stateEqLagrangian.setZero();
    stateIneqLagrangian.setZero();
    stateInputEqLagrangian.setZero();
    stateInputIneqLagrangian.setZero();
  }
};

/**
 * @brief 求一个 LagrangianMetricsGroup 的惩罚值之和。
 * @param [in] metricsGroup 拉格朗日指标组。
 * @return 所有 term 的 penalty 之和。
 */
template <typename Scalar, typename... Terms>
inline Scalar sumPenalties(
    const LagrangianMetricsGroup<Scalar, ConstraintGroupLayout<Terms...>>&
        metricsGroup) {
  return std::apply(
      [](const auto&... metrics) {
        return (Scalar(0) + ... + metrics.penalty);
      },
      metricsGroup.terms);
}

/**
 * @brief 求一个 LagrangianMetricsGroup 的约束违反平方范数之和。
 * @param [in] metricsGroup 拉格朗日指标组。
 * @return 各 term 的 constraint.squaredNorm() 之和。
 */
template <typename Scalar, typename... Terms>
inline Scalar constraintsSquaredNorm(
    const LagrangianMetricsGroup<Scalar, ConstraintGroupLayout<Terms...>>&
        metricsGroup) {
  return std::apply(
      [](const auto&... metrics) {
        return (Scalar(0) + ... + metrics.constraint.squaredNorm());
      },
      metricsGroup.terms);
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