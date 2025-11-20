#pragma once
#include "Types.hpp"
#include "ConstraintOrder.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "IntrusiveList.hpp"

/** State-input constraint function base class */
template<typename Scalar, int XDimisions, int UDimisions>
class StateInputConstraint
{
public:
  explicit StateInputConstraint(ConstraintOrder order) : order_(order) {}
  virtual ~StateInputConstraint() = default;

  /** Get the constraint order (Linear or Quadratic) */
  constexpr ConstraintOrder getOrder() const { return order_; };

  // /** Check constraint activity */
  // virtual bool isActive(scalar_t time) const { return true; }

  /** Get the constraint vector value */
  virtual Scalar getValue(const Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input) const = 0;

  /** Get the constraint linear approximation */
  virtual ScalarFunctionLinearApproximation<Scalar, XDimisions, UDimisions> getLinearApproximation(const Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input) const
  {
    static_assert(!std::is_same_v<Scalar, Scalar>, "invalid!");
  }

  /** Get the constraint quadratic approximation */
  virtual ScalarFunctionLinearApproximation<Scalar, XDimisions, UDimisions> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input) const {
    static_assert(!std::is_same_v<Scalar, Scalar>, "invalid!");
  }

private:
  ConstraintOrder order_;
};