/**
 * @file DDPSetting.hpp
 * @brief DDP/iLQR 算法参数：迭代数、收敛容差、时间步长、线搜索参数等。
 */
#pragma once
#include <cstddef>

#include "SearchStrategy/SearchStrategySettings.hpp"

/**
 * @brief DDP
 * 算法配置：最大迭代次数、代价收敛容差、时间步长及线搜索参数等。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct DDPSettings {
  /** @brief DDP 最大迭代次数。 */
  size_t maxNumIterations = 10;
  /** @brief 基于代价最小相对变化的终止条件阈值。 */
  Scalar minRelCost = Scalar(1e-3);

  /** @brief Riccati 方程积分时间步长（固定步长积分方案）。 */
  Scalar timeStep = Scalar(1e-2);

  /** @brief 线搜索策略参数。 */
  LineSearchSettings<Scalar> lineSearch{};

};  // DDP_Settings 结束。
