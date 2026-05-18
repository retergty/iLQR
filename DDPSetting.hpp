/**
 * @file DDPSetting.hpp
 * @brief DDP/iLQR 算法参数：迭代数、收敛容差、时间步长、线搜索策略等。
 */
#pragma once
#include "Integration.hpp"
#include "SearchStrategyBase.hpp"
#include "SearchStrategySettings.hpp"

/**
 * @brief DDP
 * 算法配置：最大迭代次数、代价收敛容差、约束容差、时间步长、搜索策略类型及线搜索参数等。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct DDPSettings {
  /** @brief DDP 最大迭代次数。 */
  size_t maxNumIterations_ = 10;
  /** @brief 基于代价最小相对变化的终止条件阈值。 */
  Scalar minRelCost_ = 1e-3;
  /** @brief 约束 ISE（误差平方积分）的容差。 */
  Scalar constraintTolerance_ = 1e-3;

  /** @brief Riccati 方程积分时间步长（固定步长积分方案）。 */
  Scalar timeStep_ = 1e-2;

  /** @brief 是否使用优化后的反馈策略（true）或优化后的状态-输入轨迹（false）。
   */
  bool useFeedbackPolicy_ = false;

  /** @brief 风险敏感 DDP 的风险敏感系数。 */
  Scalar riskSensitiveCoeff_ = 0.0;

  /** @brief 子问题求解策略（线搜索或 Levenberg-Marquardt）。 */
  SearchStrategyType strategy_ = SearchStrategyType::LINE_SEARCH;
  /** @brief 线搜索策略参数。 */
  LineSearchSettings<Scalar> lineSearch_{};

};  // end of DDP_Settings
