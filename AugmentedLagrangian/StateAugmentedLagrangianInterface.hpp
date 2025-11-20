#pragma once

#include "Types.hpp"
#include "Multiplier.hpp"
#include "LagrangianMetrics.hpp"
#include "QuadraticApproximation.hpp"
#include "IntrusiveList.hpp"

/** The base class for Augmented Lagrangian penalty of state constraint. */
template<typename Scalar, int XDimisions>
class StateAugmentedLagrangianInterface : IntrusiveListNode<StateAugmentedLagrangianInterface<Scalar, XDimisions>>
{
public:
  StateAugmentedLagrangianInterface() = default;
  virtual ~StateAugmentedLagrangianInterface() = default;

  /** Get the constraint and its penalty value */
  virtual LagrangianMetrics<Scalar> getValue(Scalar time, const Vector<Scalar, XDimisions>& state, const Multiplier<Scalar>& multiplier) const = 0;

  /** Get the constraint's penalty quadratic approximation */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDimisions>& state, const Multiplier<Scalar>& multiplier) const = 0;

  /** Update Lagrange/penalty multipliers and the penalty function value. */
  virtual std::pair<Multiplier<Scalar>, Scalar> updateLagrangian(const Scalar time, const Vector<Scalar, XDimisions>& state, const Scalar constraint, const Multiplier<Scalar>& multiplier) const = 0;

  /** Initialize Lagrange/penalty multipliers. */
  virtual Multiplier<Scalar> initializeLagrangian(const Scalar time) const = 0;
};