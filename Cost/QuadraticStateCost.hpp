/**
 * @file QuadraticStateCost.hpp
 * @brief 二次仅状态代价：l = 0.5 (x-x_ref)' Q (x-x_ref)，支持参考轨迹插值。
 */
#pragma once

#include "Cost.hpp"

/**
 * @brief 二次仅状态代价项：l = 0.5 (x - x_ref)' Q (x - x_ref)，x_ref 由参考轨迹插值得到。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDimisions, int ArrayLength>
class QuadraticStateCost : public StateCost<Scalar, XDimisions, ArrayLength>
{
public:
  /**
   * @brief 用权重矩阵 Q 构造二次代价。
   * @param [in] Q 半正定权重矩阵。
   */
  explicit QuadraticStateCost(const Matrix<Scalar, XDimisions, XDimisions>& Q) : StateCost<Scalar, XDimisions, ArrayLength>(0), Q_(Q) {};
  ~QuadraticStateCost() override = default;

  /** @brief 获取代价值 0.5 * (x-x_ref)' Q (x-x_ref)。 */
  Scalar getValue(Scalar time, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const final
  {
    const Vector<Scalar, XDimisions> xDeviation = getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    return 0.5 * xDeviation.dot(Q_ * xDeviation);
  }

  /** @brief 按时间索引获取代价值。 */
  Scalar getValue(
    int time_index, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const final
  {
    const Vector<Scalar, XDimisions> xDeviation = getStateDeviation(time_index, state, timeTrajectory, stateTrajectoy);
    return 0.5 * xDeviation.dot(Q_ * xDeviation);
  }

  /** @brief 获取代价的二次近似（dfdxx=Q, dfdx=Q*(x-x_ref), f=0.5*(x-x_ref)'*dfdx）。 */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
    getQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDimisions>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const final
  {
    const Vector<Scalar, XDimisions> xDeviation = getStateDeviation(time, state, timeTrajectory, stateTrajectoy);

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
    const Vector<Scalar, XDimisions> xDeviation = getStateDeviation(time_index, state, timeTrajectory, stateTrajectoy);

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
      (void)timeTrajectory;
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
    const Matrix<Scalar, UDimisions, XDimisions>& P) : StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength>(0), Q_(Q), R_(R), P_(P)
  {
    has_P_ = true;
  };

  QuadraticStateInputCost(const Matrix<Scalar, XDimisions, XDimisions>& Q,
    const Matrix<Scalar, UDimisions, UDimisions>& R) : StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength>(0), Q_(Q), R_(R)
  {
    P_.setZero();
    has_P_ = false;
  };

  ~QuadraticStateInputCost() override = default;

  /** Get cost term value */
  Scalar getValue(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const final
  {
    Vector<Scalar, XDimisions> stateDeviation = getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    Vector<Scalar, UDimisions> inputDeviation = getInputDeviation(time, input, timeTrajectory, inputTrajectory);

    if (has_P_)
    {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) + 0.5 * inputDeviation.dot(R_ * inputDeviation) +
      inputDeviation.dot(P_ * stateDeviation);
    }
    else
    {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) + 0.5 * inputDeviation.dot(R_ * inputDeviation);
    }
  }

  /** Get cost term value */
  Scalar getValue(int time_index, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const final
  {
    (void)timeTrajectory;
    Vector<Scalar, XDimisions> stateDeviation = getStateDeviation(time_index, state, timeTrajectory, stateTrajectoy);
    Vector<Scalar, UDimisions> inputDeviation = getInputDeviation(time_index, input, timeTrajectory, inputTrajectory);

    if (has_P_)
    {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) + 0.5 * inputDeviation.dot(R_ * inputDeviation) +
        inputDeviation.dot(P_ * stateDeviation);
    }
    else
    {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) + 0.5 * inputDeviation.dot(R_ * inputDeviation);
    }
  }

  /** Get cost term quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>
    getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const final
  {
      (void)timeTrajectory;
    Vector<Scalar, XDimisions> stateDeviation = getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    Vector<Scalar, UDimisions> inputDeviation = getInputDeviation(time, input, timeTrajectory, inputTrajectory);

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
    return state - LinearInterpolation::interpolate(time, timeTrajectory, stateTrajectoy);
  }
  /** Computes the state-input deviation pair around the nominal state and input.
   * This method can be overwritten if desiredTrajectory has a different dimensions. */
  Vector<Scalar, XDimisions> getStateDeviation(int time_index, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const
  {
    (void)timeTrajectory;
    return  state - stateTrajectoy[time_index];
  }
  Vector<Scalar, UDimisions> getInputDeviation(Scalar time, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const
  {
    return input - LinearInterpolation::interpolate(time, timeTrajectory, inputTrajectory);
  }
  Vector<Scalar, UDimisions> getInputDeviation(int time_index, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const
  {
    (void)timeTrajectory;
    return input - inputTrajectory[time_index];
  }

private:
  Matrix<Scalar, XDimisions, XDimisions> Q_;
  Matrix<Scalar, UDimisions, UDimisions> R_;
  Matrix<Scalar, UDimisions, XDimisions> P_;
  bool has_P_{ true };
};
