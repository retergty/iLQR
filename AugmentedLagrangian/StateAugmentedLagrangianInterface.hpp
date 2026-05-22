/**
 * @file StateAugmentedLagrangianInterface.hpp
 * @brief 仅状态增广拉格朗日接口：约束惩罚取值、二次近似与乘子初始化/更新。
 */
#pragma once

#include "LagrangianMetrics.hpp"
#include "Multiplier.hpp"
#include "QuadraticApproximation.hpp"
#include "Types.hpp"

/**
 * @brief
 * 仅状态约束的增广拉格朗日惩罚接口：提供取值、二次近似、乘子更新与初始化。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 */
template <typename Scalar, int XDim, int CDim>
class StateAugmentedLagrangianInterface {
 public:
  StateAugmentedLagrangianInterface() = default;
  virtual ~StateAugmentedLagrangianInterface() = default;

  /** @brief 获取约束与惩罚值（LagrangianMetrics）。 */
  virtual LagrangianMetrics<Scalar, CDim> getValue(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Multiplier<Scalar, CDim>& multiplier) const = 0;

  /** @brief 获取约束惩罚的二次近似。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>
  getQuadraticApproximation(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Multiplier<Scalar, CDim>& multiplier) const = 0;

  /** @brief 更新拉格朗日/惩罚乘子并返回更新后乘子与惩罚值。 */
  virtual std::pair<Multiplier<Scalar, CDim>, Scalar> updateLagrangian(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, CDim>& constraint,
      const Multiplier<Scalar, CDim>& multiplier) const = 0;

  /** @brief 初始化拉格朗日/惩罚乘子。 */
  virtual Multiplier<Scalar, CDim> initializeLagrangian(
      const Scalar time) const = 0;
};