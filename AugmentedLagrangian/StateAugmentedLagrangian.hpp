#pragma once
#include "Types.hpp"
#include "StateAugmentedLagrangianInterface.hpp"
#include "StateConstraint.hpp"
#include "AugmentedPenaltyBase.hpp"
#include "Penalty.hpp"
#include "LagrangianMetrics.hpp"
#include "QuadraticApproximation.hpp"

/** The base class for Augmented Lagrangian penalty of state constraint. */
template <typename Scalar, int XDimisions>
class StateAugmentedLagrangian final : public StateAugmentedLagrangianInterface<Scalar, XDimisions>
{
public:
  /**
   * Constructor.
   * @param [in] constraintPtr: A pointer to the constraint which will be enforced as soft constraints.
   * @param [in] augmented_penalty: pointer to the penalty function on the constraint.
   */
  StateAugmentedLagrangian(StateConstraint<Scalar, XDimisions> *constraint, AugmentedPenaltyBase<Scalar> *augmented_penalty) : constraint_ptr_(constraint), penalty_(augmented_penalty) {};

  LagrangianMetrics<Scalar> getValue(const Scalar time, const Vector<Scalar, XDimisions> &state, const Multiplier<Scalar> &multiplier) const override
  {
    const Scalar h = constraint_ptr_->getValue(time, state);
    const Scalar p = multiplier.penalty * penalty_.getValue(time, h, &multiplier.lagrangian);
    return {p, h};
  }

  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDimisions> &state, const Multiplier<Scalar> &multiplier) const override
  {
    switch (constraint_ptr_->getOrder())
    {
    case ConstraintOrder::Linear:
      return multiplier.penalty * penalty_.getQuadraticApproximation(time, constraint_ptr_->getLinearApproximation(time, state), &multiplier.lagrangian);
    case ConstraintOrder::Quadratic:
      return multiplier.penalty * penalty_.getQuadraticApproximation(time, constraint_ptr_->getQuadraticApproximation(time, state), &multiplier.lagrangian);
    default:
      return ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>();
    }
  }

  std::pair<Multiplier<Scalar>, Scalar> updateLagrangian(const Scalar time, const Vector<Scalar, XDimisions> &state, const Scalar constraint, const Multiplier<Scalar> &multiplier) const override
  {
    (void)state;
    const Multiplier<Scalar> updatedMultiplier{multiplier.penalty, penalty_.updateMultipliers(time, constraint, multiplier.lagrangian)};
    const Scalar penalty = updatedMultiplier.penalty * penalty_.getValue(time, constraint, updatedMultiplier.lagrangian);
    return {updatedMultiplier, penalty};
  }
  Multiplier<Scalar> initializeLagrangian(const Scalar time) const override
  {
    return {1.0, penalty_.initializeMultipliers()};
  }

private:
  StateConstraint<Scalar, XDimisions> *constraint_ptr_;
  Penalty<Scalar, XDimisions, 0> penalty_;
};