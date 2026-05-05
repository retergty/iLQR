/**
 * @file StateAugmentedLagrangian.hpp
 * @brief 仅状态增广拉格朗日实现：绑定约束与惩罚，提供取值、二次近似与乘子更新。
 */
#pragma once
#include "Types.hpp"
#include "StateAugmentedLagrangianInterface.hpp"
#include "StateConstraint.hpp"
#include "AugmentedPenaltyBase.hpp"
#include "Penalty.hpp"
#include "LagrangianMetrics.hpp"
#include "QuadraticApproximation.hpp"

/**
 * @brief 仅状态约束的增广拉格朗日惩罚实现：委托约束与 Penalty 计算取值与二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 */
template <typename Scalar, int XDim>
class StateAugmentedLagrangian final : public StateAugmentedLagrangianInterface<Scalar, XDim>
{
public:
  /**
   * @brief 用约束指针与惩罚指针构造。
   * @param [in] constraint 作为软约束的仅状态约束。
   * @param [in] augmented_penalty 约束上的惩罚函数。
   */
  StateAugmentedLagrangian(StateConstraint<Scalar, XDim> *constraint, AugmentedPenaltyBase<Scalar> *augmented_penalty) : constraint_ptr_(constraint), penalty_(augmented_penalty) {};

  LagrangianMetrics<Scalar> getValue(const Scalar time, const Vector<Scalar, XDim> &state, const Multiplier<Scalar> &multiplier) const override
  {
    const Scalar h = constraint_ptr_->getValue(time, state);
    const Scalar p = multiplier.penalty * penalty_.getValue(time, h, &multiplier.lagrangian);
    return {p, h};
  }

  ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDim> &state, const Multiplier<Scalar> &multiplier) const override
  {
    switch (constraint_ptr_->getOrder())
    {
    case ConstraintOrder::Linear:
      return multiplier.penalty * penalty_.getQuadraticApproximation(time, constraint_ptr_->getLinearApproximation(time, state), &multiplier.lagrangian);
    case ConstraintOrder::Quadratic:
      return multiplier.penalty * penalty_.getQuadraticApproximation(time, constraint_ptr_->getQuadraticApproximation(time, state), &multiplier.lagrangian);
    default:
      return ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>();
    }
  }

  std::pair<Multiplier<Scalar>, Scalar> updateLagrangian(const Scalar time, const Vector<Scalar, XDim> &state, const Scalar constraint, const Multiplier<Scalar> &multiplier) const override
  {
    (void)state;
    const Multiplier<Scalar> updatedMultiplier{multiplier.penalty, penalty_.updateMultipliers(time, constraint, multiplier.lagrangian)};
    const Scalar penalty = updatedMultiplier.penalty * penalty_.getValue(time, constraint, updatedMultiplier.lagrangian);
    return {updatedMultiplier, penalty};
  }
  Multiplier<Scalar> initializeLagrangian(const Scalar time) const override
  {
    (void)time;
    return {1.0, penalty_.initializeMultipliers()};
  }

private:
  StateConstraint<Scalar, XDim> *constraint_ptr_;
  Penalty<Scalar, XDim, 0> penalty_;
};