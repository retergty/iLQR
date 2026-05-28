/**
 * @file StateInputAugmentedLagrangianInterface.hpp
 * @brief 状态-输入增广拉格朗日接口：约束惩罚取值、二次近似与乘子初始化/更新。
 */
#pragma once

#include "Approximation/QuadraticApproximation.hpp"
#include "Matrix/Types.hpp"
#include "ModelData/Multiplier.hpp"
#include "OptimalControlData/LagrangianMetrics.hpp"

/**
 * @brief
 * 状态-输入约束的增广拉格朗日惩罚接口：提供取值、二次近似、乘子更新与初始化。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int XDim, int UDim, int CDim>
class StateInputAugmentedLagrangianInterface {
 public:
  StateInputAugmentedLagrangianInterface() = default;
  virtual ~StateInputAugmentedLagrangianInterface() = default;

  /** 获取约束及其惩罚值。 */
  virtual LagrangianMetrics<Scalar, CDim> getValue(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const Multiplier<Scalar, CDim>& lagrangian) const = 0;

  /** 获取约束惩罚的二次近似。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const Multiplier<Scalar, CDim>& lagrangian) const = 0;

  /** 更新拉格朗日/惩罚乘子以及惩罚函数值。 */
  virtual std::pair<Multiplier<Scalar, CDim>, Scalar> updateLagrangian(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input, const Vector<Scalar, CDim>& constraint,
      const Multiplier<Scalar, CDim>& lagrangian) const = 0;

  /** 初始化拉格朗日/惩罚乘子。 */
  virtual Multiplier<Scalar, CDim> initializeLagrangian(Scalar time) const = 0;
};
