/**
 * @file SearchStrategySettings.hpp
 * @brief 搜索策略参数：线搜索设置。
 */
#pragma once
#include "iLQR/HessianCorrection.hpp"

/**
 * @brief 搜索策略基础设置：代价相对变化收敛阈值。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct SearchStrategyBaseSettings {
  /** @brief 基于代价最小相对变化的终止条件阈值。 */
  Scalar minRelCost{Scalar(1e-3)};
  /** @brief 计算代价相对变化时使用的最小归一化参考量。 */
  Scalar costNormalizationBase{Scalar(1.0)};
};

/**
 * @brief 线搜索策略参数：步长范围、收缩率、Armijo 系数与 Hessian 修正。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct LineSearchSettings : public SearchStrategyBaseSettings<Scalar> {
  /** @brief 线搜索最小步长。 */
  Scalar minStepLength = Scalar(0.05);
  /** @brief 线搜索最大步长。 */
  Scalar maxStepLength = Scalar(1.0);
  /** @brief 步长收缩率。 */
  Scalar contractionRate = Scalar(0.5);
  /** @brief Armijo 系数：满足 merit(u+a*p) < merit(u) + c*a*dfdu.dot(p)
   * 时接受。 */
  Scalar armijoCoefficient = Scalar(1e-4);
  /** @brief Hessian 修正策略（如对角平移）。 */
  HessianCorrectionStrategy hessianCorrectionStrategy =
      HessianCorrectionStrategy::DIAGONAL_SHIFT;
  /** @brief Riccati 反向递推数值稳定性用的 Hessian 修正倍数。 */
  Scalar hessianCorrectionMultiple = Scalar(1e-6);
};
