/**
 * @file RolloutTraits.hpp
 * @brief 根据动力学模式在编译期选择连续或离散 rollout 实现。
 */
#pragma once

#include "Rollout/DiscreteTimeRollout.hpp"
#include "Rollout/TimeTriggeredRollout.hpp"
#include "iLQR/iLQRDescriptor.hpp"

/**
 * @brief rollout 模式 traits。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @tparam DynamicsMode 动力学模式 tag。
 */
template <typename Scalar, int XDim, int UDim, typename DynamicsMode>
struct RolloutModeTraits;

/** @brief 连续时间动力学使用 RK 积分 rollout。 */
template <typename Scalar, int XDim, int UDim>
struct RolloutModeTraits<Scalar, XDim, UDim, ContinuousDynamics> {
  using Rollout_t = TimeTriggeredRollout<Scalar, XDim, UDim>;
};

/** @brief 离散时间动力学使用直接状态转移 rollout。 */
template <typename Scalar, int XDim, int UDim>
struct RolloutModeTraits<Scalar, XDim, UDim, DiscreteDynamics> {
  using Rollout_t = DiscreteTimeRollout<Scalar, XDim, UDim>;
};
