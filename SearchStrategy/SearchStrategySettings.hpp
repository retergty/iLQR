/**
 * @file SearchStrategySettings.hpp
 * @brief 搜索策略参数：DDP 策略类型、线搜索与 Levenberg-Marquardt 设置。
 */
#pragma once
#include "iLQR/HessianCorrection.hpp"

/**
 * @brief DDP 子问题求解策略枚举：线搜索或 Levenberg-Marquardt。
 */
enum class SearchStrategyType { LINE_SEARCH, LEVENBERG_MARQUARDT };

/**
 * @brief 搜索策略基础设置：日志、调试输出与代价相对变化收敛阈值。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct SearchStrategyBaseSettings {
  /** @brief 是否输出 DDP 日志。 */
  bool displayInfo = false;
  /** @brief 是否打印 rollout 轨迹用于调试。 */
  bool debugPrintRollout = false;
  /** @brief 基于代价最小相对变化的终止条件阈值。 */
  Scalar minRelCost{1e-3};
};

/**
 * @brief 线搜索策略参数：步长范围、收缩率、Armijo 系数与 Hessian 修正。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct LineSearchSettings : public SearchStrategyBaseSettings<Scalar> {
  /** @brief 线搜索最小步长。 */
  Scalar minStepLength = 0.05;
  /** @brief 线搜索最大步长。 */
  Scalar maxStepLength = 1.0;
  /** @brief 步长收缩率。 */
  Scalar contractionRate = 0.5;
  /** @brief Armijo 系数：满足 merit(u+a*p) < merit(u) + c*a*dfdu.dot(p)
   * 时接受。 */
  Scalar armijoCoefficient = 1e-4;
  /** @brief Hessian 修正策略（如对角平移）。 */
  HessianCorrectionStrategy hessianCorrectionStrategy =
      HessianCorrectionStrategy::DIAGONAL_SHIFT;
  /** @brief Riccati 反向递推数值稳定性用的 Hessian 修正倍数。 */
  Scalar hessianCorrectionMultiple = 1e-6;
};

/**
 * @brief Levenberg-Marquardt 策略参数：接受比、Riccati
 * 倍数递推与连续拒绝次数上限。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct LevenbergMarquardtSettings : public SearchStrategyBaseSettings<Scalar> {
  /** @brief 接受迭代解的最小 rho（实际下降/预测下降），应在 [0, 0.25)。 */
  Scalar minAcceptedPho = 0.25;
  /** @brief Riccati 倍数几何递推的默认比例。 */
  Scalar riccatiMultipleDefaultRatio = 2.0;
  /** @brief Riccati 倍数几何递推的默认因子。 */
  Scalar riccatiMultipleDefaultFactor = 1e-6;
  /** @brief 连续拒绝迭代解的最大次数。 */
  int maxNumSuccessiveRejections = 5;
};