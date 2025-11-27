#pragma once
#include "Types.hpp"
#include "Reference.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "IntrusiveList.hpp"
#include "LinearInterpolation.hpp"
/** State-only cost term */
template <typename Scalar, int XDimisions, int ArrayLength>
class StateCost : public IntrusiveListNode<StateCost<Scalar, XDimisions, ArrayLength>>
{
public:
  StateCost(int cost_number) : number(cost_number) {};
  virtual ~StateCost() = default;

  /** Check if cost term is active */
  virtual bool isActive(Scalar time) const { return true; }

  /** Get cost term value */
  virtual Scalar getValue(Scalar time, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const = 0;

  /** Get cost term quadratic approximation */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
    getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const = 0;

  // indentify state cost, must be unique
  int number;

protected:
  StateCost(const StateCost& rhs) = default;
};

/** State-input cost term */
template <typename Scalar, int XDimisions, int UDimisions, int ArrayLength>
class StateInputCost : public IntrusiveListNode<StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength>>
{
public:
  StateInputCost(int cost_number) : number(cost_number) {};
  virtual ~StateInputCost() = default;

  /** Check if cost term is active */
  virtual bool isActive(Scalar time) const { return true; }

  /** Get cost term value */
  virtual Scalar getValue(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const = 0;

  /** Get cost term quadratic approximation */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>
    getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const = 0;

  // identify state input cost, must be unique
  int number;

protected:
  StateInputCost(const StateInputCost& rhs) = default;
};