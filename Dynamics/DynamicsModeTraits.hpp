/**
 * @file DynamicsModeTraits.hpp
 * @brief 根据动力学模式在编译期选择连续或离散动力学基类。
 */
#pragma once

#include "DiscreteSystemDynamicsBase.hpp"
#include "SystemDynamicsBase.hpp"
#include "iLQRDescriptor.hpp"

/**
 * @brief 动力学模式 traits。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @tparam DynamicsMode 动力学模式 tag。
 */
template <typename Scalar, int XDim, int UDim, typename DynamicsMode>
struct DynamicsModeTraits;

/** @brief 连续时间动力学使用 SystemDynamicsBase。 */
template <typename Scalar, int XDim, int UDim>
struct DynamicsModeTraits<Scalar, XDim, UDim, ContinuousDynamics> {
  using DynamicsBase_t = SystemDynamicsBase<Scalar, XDim, UDim>;
};

/** @brief 离散时间动力学使用 DiscreteSystemDynamicsBase。 */
template <typename Scalar, int XDim, int UDim>
struct DynamicsModeTraits<Scalar, XDim, UDim, DiscreteDynamics> {
  using DynamicsBase_t = DiscreteSystemDynamicsBase<Scalar, XDim, UDim>;
};
