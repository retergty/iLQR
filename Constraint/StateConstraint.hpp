#pragma once

#include "Types.hpp"
#include "ConstraintOrder.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include <type_traits>

/** State-only constraint function base class */
template<typename Scalar, int XDimisions>
class StateConstraint
{
public:
  explicit StateConstraint(ConstraintOrder order) : order_(order) {}
  virtual ~StateConstraint() = default;

  /** Get the constraint order (Linear or Quadratic) */
  constexpr ConstraintOrder getOrder() const { return order_; };

  // /** Check constraint activity */
  // virtual bool isActive(const Scalar time) const { return true; }

  /** Get the constraint value */
  virtual Scalar getValue(const Scalar time, const Vector<Scalar, XDimisions>& state) const = 0;

  /** Get the constraint linear approximation */
  virtual ScalarFunctionLinearApproximation<Scalar, XDimisions, 0> getLinearApproximation(const Scalar time, const Vector<Scalar, XDimisions>& state) const
  {
    static_assert(!std::is_same_v<Scalar, Scalar>, "invalid!");
  }

  /** Get the constraint quadratic approximation */
  virtual ScalarFunctionLinearApproximation<Scalar, XDimisions, 0> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDimisions>& state) const
  {
    static_assert(!std::is_same_v<Scalar, Scalar>, "invalid!");
  }

private:
  ConstraintOrder order_;
};
