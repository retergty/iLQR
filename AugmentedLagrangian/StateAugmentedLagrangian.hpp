/**
 * @file StateAugmentedLagrangian.hpp
 * @brief 仅状态增广拉格朗日实现：绑定约束与惩罚，提供取值、二次近似与乘子更新。
 */
#pragma once
#include "AugmentedPenaltyBase.hpp"
#include "LagrangianMetrics.hpp"
#include "MultidimensionalPenalty.hpp"
#include "QuadraticApproximation.hpp"
#include "StateAugmentedLagrangianInterface.hpp"
#include "StateConstraint.hpp"
#include "Types.hpp"

/**
 * @brief 仅状态约束的增广拉格朗日惩罚实现：委托约束与 Penalty
 * 计算取值与二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 */
template <typename Scalar, int XDim, int CDim>
class StateAugmentedLagrangian final
    : public StateAugmentedLagrangianInterface<Scalar, XDim, CDim> {
 public:
  /**
   * @brief 用约束指针与惩罚指针构造。
   * @param [in] constraint 作为软约束的仅状态约束。
   * @param [in] augmented_penalty 约束上的惩罚函数。
   */
  StateAugmentedLagrangian(StateConstraint<Scalar, XDim, CDim>* constraint,
                           AugmentedPenaltyBase<Scalar>* augmented_penalty)
      : constraintPtr_(constraint), penalty_(augmented_penalty) {};

  LagrangianMetrics<Scalar, CDim> getValue(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Multiplier<Scalar, CDim>& multiplier) const override {
    const Vector<Scalar, CDim> h = constraintPtr_->getValue(time, state);
    const Scalar p =
        multiplier.penalty * penalty_.getValue(time, h, multiplier.lagrangian);
    return {p, h};
  }

  ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>
  getQuadraticApproximation(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Multiplier<Scalar, CDim>& multiplier) const override {
    switch (constraintPtr_->getOrder()) {
      case ConstraintOrder::Linear:
        return multiplier.penalty *
               penalty_.getQuadraticApproximation(
                   time, constraintPtr_->getLinearApproximation(time, state),
                   multiplier.lagrangian);
      case ConstraintOrder::Quadratic:
        return multiplier.penalty *
               penalty_.getQuadraticApproximation(
                   time, constraintPtr_->getQuadraticApproximation(time, state),
                   multiplier.lagrangian);
      default:
        return ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>::Zero();
    }
  }

  std::pair<Multiplier<Scalar, CDim>, Scalar> updateLagrangian(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, CDim>& constraint,
      const Multiplier<Scalar, CDim>& multiplier) const override {
    (void)state;
    const Multiplier<Scalar, CDim> updatedMultiplier{
        multiplier.penalty,
        penalty_.updateMultipliers(time, constraint, multiplier.lagrangian)};
    const Scalar penalty =
        updatedMultiplier.penalty *
        penalty_.getValue(time, constraint, updatedMultiplier.lagrangian);
    return {updatedMultiplier, penalty};
  }
  Multiplier<Scalar, CDim> initializeLagrangian(
      const Scalar time) const override {
    (void)time;
    return {static_cast<Scalar>(1), penalty_.initializeMultipliers()};
  }

 private:
  StateConstraint<Scalar, XDim, CDim>* constraintPtr_;
  MultidimensionalPenalty<Scalar, XDim, 0, CDim> penalty_;
};