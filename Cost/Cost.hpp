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

template <typename Scalar, int XDimisions, int ArrayLength>
/** Quadratic state-only cost term */
class QuadraticStateCost : public StateCost<Scalar, XDimisions, ArrayLength>
{
public:
  /**
   * Constructor for the quadratic cost function defined as the following:
   * \f$ \l = 0.5(x-x_{n})' Q (x-x_{n}) \f$.
   * @param [in] Q: \f$ Q \f$
   */
  explicit QuadraticStateCost(const Matrix<Scalar, XDimisions, XDimisions>& Q) : Q_(Q) {};
  ~QuadraticStateCost() override = default;

  /** Get cost term value */
  Scalar getValue(Scalar time, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const final
  {
    const Vector<Scalar, XDimisions> xDeviation = getStateDeviation(time, state, stateTrajectoy);
    return 0.5 * xDeviation.dot(Q_ * xDeviation);
  }

  /** Get cost term value */
  Scalar getValue(
    int time_index, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const final
  {
    const Vector<Scalar, XDimisions> xDeviation = getStateDeviation(time_index, state, stateTrajectoy);
    return 0.5 * xDeviation.dot(Q_ * xDeviation);
  }

  /** Get cost term quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
    getQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDimisions>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const final
  {
    const Vector<Scalar, XDimisions> xDeviation = getStateDeviation(time, state, stateTrajectoy);

    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> Phi;
    Phi.dfdxx = Q_;
    Phi.dfdx = Q_ * xDeviation;
    Phi.f = 0.5 * xDeviation.dot(Phi.dfdx);
    return Phi;
  }

  /** Get cost term quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
    getQuadraticApproximation(
      int time_index, const Vector<Scalar, XDimisions>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const final
  {
    const Vector<Scalar, XDimisions> xDeviation = getStateDeviation(time_index, state, stateTrajectoy);

    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> Phi;
    Phi.dfdxx = Q_;
    Phi.dfdx = Q_ * xDeviation;
    Phi.f = 0.5 * xDeviation.dot(Phi.dfdx);
    return Phi;
  }

protected:
  QuadraticStateCost(const QuadraticStateCost& rhs) = default;

  /** Computes the state deviation for the nominal state.
   * This method can be overwritten if desiredTrajectory has a different dimensions. */
  Vector<Scalar, XDimisions> getStateDeviation(Scalar time, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const
  {

    return state - LinearInterpolation::interpolate(time, timeTrajectory, stateTrajectoy);
  }

  /** Computes the state deviation for the nominal state.
   * This method can be overwritten if desiredTrajectory has a different dimensions. */
  Vector<Scalar, XDimisions> getStateDeviation(int time_index, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const
  {
    return state - stateTrajectoy[time_index];
  }

private:
  Matrix<Scalar, XDimisions, XDimisions> Q_;
};

/** Quadratic state-input cost term */
template <typename Scalar, int XDimisions, int UDimisions, int ArrayLength>
class QuadraticStateInputCost : public StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength>
{
public:
  /**
   * Constructor for the quadratic cost function defined as the following:
   * \f$ L = 0.5(x-x_{n})' Q (x-x_{n}) + 0.5(u-u_{n})' R (u-u_{n}) + (u-u_{n})' P (x-x_{n}) \f$
   * @param [in] Q: \f$ Q \f$
   * @param [in] R: \f$ R \f$
   * @param [in] P: \f$ P \f$
   */
  QuadraticStateInputCost(const Matrix<Scalar, XDimisions, XDimisions>& Q,
    const Matrix<Scalar, UDimisions, UDimisions>& R,
    const Matrix<Scalar, UDimisions, XDimisions>& P) : Q_(Q), R_(R), P_(P)
  {
    has_P_ = true;
  };

  QuadraticStateInputCost(const Matrix<Scalar, XDimisions, XDimisions>& Q,
    const Matrix<Scalar, UDimisions, UDimisions>& R)
  {
    P_.setZero();
    has_P_ = false;
  };

  ~QuadraticStateInputCost() override = default;

  /** Get cost term value */
  Scalar getValue(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const final
  {
    Vector<Scalar, XDimisions> stateDeviation = getStateDeviation(time, state, stateTrajectoy);
    Vector<Scalar, UDimisions> inputDeviation = getInputDeviation(time, input, inputTrajectory);

    if (has_P_)
    {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) + 0.5 * inputDeviation.dot(R_ * inputDeviation);
    }
    else
    {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) + 0.5 * inputDeviation.dot(R_ * inputDeviation) +
        inputDeviation.dot(P_ * stateDeviation);
    }
  }

  /** Get cost term value */
  Scalar getValue(int time_index, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const final
  {
    Vector<Scalar, XDimisions> stateDeviation = getStateDeviation(time_index, state, stateTrajectoy);
    Vector<Scalar, UDimisions> inputDeviation = getInputDeviation(time_index, input, inputTrajectory);

    if (has_P_)
    {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) + 0.5 * inputDeviation.dot(R_ * inputDeviation);
    }
    else
    {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) + 0.5 * inputDeviation.dot(R_ * inputDeviation) +
        inputDeviation.dot(P_ * stateDeviation);
    }
  }

  /** Get cost term quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>
    getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const final
  {
    Vector<Scalar, XDimisions> stateDeviation = getStateDeviation(time, state, stateTrajectoy);
    Vector<Scalar, UDimisions> inputDeviation = getInputDeviation(time, input, inputTrajectory);

    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> L;
    L.dfdxx = Q_;
    L.dfduu = R_;
    L.dfdx = Q_ * stateDeviation;
    L.dfdu = R_ * inputDeviation;
    L.f = 0.5 * stateDeviation.dot(L.dfdx) + 0.5 * inputDeviation.dot(L.dfdu);

    if (has_P_ == 0)
    {
      L.dfdux.setZero();
    }
    else
    {
      const Vector<Scalar, UDimisions> pDeviation = P_ * stateDeviation;
      L.f += inputDeviation.dot(pDeviation);
      L.dfdu += pDeviation;
      L.dfdx += P_.transpose() * inputDeviation;
      L.dfdux = P_;
    }

    return L;
  }

protected:
  QuadraticStateInputCost(const QuadraticStateInputCost& rhs) = default;

  /** Computes the state-input deviation pair around the nominal state and input.
   * This method can be overwritten if desiredTrajectory has a different dimensions. */
  Vector<Scalar, XDimisions> getStateDeviation(Scalar time, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const
  {
    return state - LinearInterpolation::interpolate(time, timeTrajectory, timeTrajectory);
  }
  /** Computes the state-input deviation pair around the nominal state and input.
   * This method can be overwritten if desiredTrajectory has a different dimensions. */
  Vector<Scalar, XDimisions> getStateDeviation(int time_index, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const
  {
    return  state - stateTrajectoy[time_index];
  }
  Vector<Scalar, XDimisions> getInputDeviation(Scalar time, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const
  {
    return input - LinearInterpolation::interpolate(time, timeTrajectory, inputTrajectory);
  }
  Vector<Scalar, XDimisions> getInputDeviation(int time_index, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const
  {
    return input - inputTrajectory[time_index];
  }

private:
  Matrix<Scalar, XDimisions, XDimisions> Q_;
  Matrix<Scalar, UDimisions, UDimisions> R_;
  Matrix<Scalar, UDimisions, XDimisions> P_;
  bool has_P_{ true };
};
