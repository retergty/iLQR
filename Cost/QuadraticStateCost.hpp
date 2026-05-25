/**
 * @file QuadraticStateCost.hpp
 * @brief 二次仅状态代价：l = 0.5 (x-x_ref)' Q (x-x_ref)，支持参考轨迹插值。
 */
#pragma once

#include "Cost.hpp"
#include "LinearInterpolation.hpp"

/**
 * @brief 二次仅状态代价项：l = 0.5 (x - x_ref)' Q (x - x_ref)，x_ref
 * 由参考轨迹插值得到。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDim, int ArrayLength>
class QuadraticStateCost : public StateCost<Scalar, XDim, ArrayLength> {
 public:
  /**
   * @brief 用权重矩阵 Q 构造二次代价。
   * @param [in] Q 半正定权重矩阵。
   */
  explicit QuadraticStateCost(const Matrix<Scalar, XDim, XDim>& Q)
      : StateCost<Scalar, XDim, ArrayLength>(0), Q_(Q) {};
  ~QuadraticStateCost() override = default;

  /** @brief 获取代价值 0.5 * (x-x_ref)' Q (x-x_ref)。 */
  Scalar getValue(Scalar time, const Vector<Scalar, XDim>& state,
                  const std::array<Scalar, ArrayLength>& timeTrajectory,
                  const std::array<Vector<Scalar, XDim>, ArrayLength>&
                      stateTrajectoy) const final {
    const Vector<Scalar, XDim> xDeviation =
        getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    return 0.5 * xDeviation.dot(Q_ * xDeviation);
  }

  /** @brief 获取代价的二次近似（dfdxx=Q, dfdx=Q*(x-x_ref),
   * f=0.5*(x-x_ref)'*dfdx）。 */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy)
      const final {
    const Vector<Scalar, XDim> xDeviation =
        getStateDeviation(time, state, timeTrajectory, stateTrajectoy);

    ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> Phi;
    Phi.dfdxx = Q_;
    Phi.dfdx = Q_ * xDeviation;
    Phi.f = 0.5 * xDeviation.dot(Phi.dfdx);
    return Phi;
  }

 protected:
  QuadraticStateCost(const QuadraticStateCost& rhs) = default;

  /** 计算相对于名义状态的状态偏差。
   * 如果 desiredTrajectory 具有不同
   * 维度，可重写此方法。 */
  Vector<Scalar, XDim> getStateDeviation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy)
      const {
    return state - LinearInterpolation::interpolate(time, timeTrajectory,
                                                    stateTrajectoy);
  }

 private:
  Matrix<Scalar, XDim, XDim> Q_;
};

/** 二次状态-输入代价项。 */
template <typename Scalar, int XDim, int UDim, int ArrayLength>
class QuadraticStateInputCost
    : public StateInputCost<Scalar, XDim, UDim, ArrayLength> {
 public:
  /**
   * 构造以下形式定义的二次代价函数：
   * \f$ L = 0.5(x-x_{n})' Q (x-x_{n}) + 0.5(u-u_{n})' R (u-u_{n}) + (u-u_{n})'
   * P (x-x_{n}) \f$
   * @param [in] Q: \f$ Q \f$
   * @param [in] R: \f$ R \f$
   * @param [in] P: \f$ P \f$
   */
  QuadraticStateInputCost(const Matrix<Scalar, XDim, XDim>& Q,
                          const Matrix<Scalar, UDim, UDim>& R,
                          const Matrix<Scalar, UDim, XDim>& P)
      : StateInputCost<Scalar, XDim, UDim, ArrayLength>(0),
        Q_(Q),
        R_(R),
        P_(P) {
    has_P_ = true;
  };

  QuadraticStateInputCost(const Matrix<Scalar, XDim, XDim>& Q,
                          const Matrix<Scalar, UDim, UDim>& R)
      : StateInputCost<Scalar, XDim, UDim, ArrayLength>(0), Q_(Q), R_(R) {
    P_.setZero();
    has_P_ = false;
  };

  ~QuadraticStateInputCost() override = default;

  /** 获取代价项取值。 */
  Scalar getValue(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory)
      const final {
    Vector<Scalar, XDim> stateDeviation =
        getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    Vector<Scalar, UDim> inputDeviation =
        getInputDeviation(time, input, timeTrajectory, inputTrajectory);

    if (has_P_) {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) +
             0.5 * inputDeviation.dot(R_ * inputDeviation) +
             inputDeviation.dot(P_ * stateDeviation);
    } else {
      return 0.5 * stateDeviation.dot(Q_ * stateDeviation) +
             0.5 * inputDeviation.dot(R_ * inputDeviation);
    }
  }

  /** 获取代价项二次近似。 */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory)
      const final {
    (void)timeTrajectory;
    Vector<Scalar, XDim> stateDeviation =
        getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    Vector<Scalar, UDim> inputDeviation =
        getInputDeviation(time, input, timeTrajectory, inputTrajectory);

    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> L;
    L.dfdxx = Q_;
    L.dfduu = R_;
    L.dfdx = Q_ * stateDeviation;
    L.dfdu = R_ * inputDeviation;
    L.f = 0.5 * stateDeviation.dot(L.dfdx) + 0.5 * inputDeviation.dot(L.dfdu);

    if (has_P_ == 0) {
      L.dfdux.setZero();
    } else {
      const Vector<Scalar, UDim> pDeviation = P_ * stateDeviation;
      L.f += inputDeviation.dot(pDeviation);
      L.dfdu += pDeviation;
      L.dfdx += P_.transpose() * inputDeviation;
      L.dfdux = P_;
    }

    return L;
  }

 protected:
  QuadraticStateInputCost(const QuadraticStateInputCost& rhs) = default;

  /** 计算相对于名义状态和输入的状态-输入偏差对。
   * 如果 desiredTrajectory 具有不同维度，可重写此方法。
   * 维度，可重写此方法。 */
  Vector<Scalar, XDim> getStateDeviation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy)
      const {
    return state - LinearInterpolation::interpolate(time, timeTrajectory,
                                                    stateTrajectoy);
  }
  Vector<Scalar, UDim> getInputDeviation(
      Scalar time, const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory)
      const {
    return input - LinearInterpolation::interpolate(time, timeTrajectory,
                                                    inputTrajectory);
  }

 private:
  Matrix<Scalar, XDim, XDim> Q_;
  Matrix<Scalar, UDim, UDim> R_;
  Matrix<Scalar, UDim, XDim> P_;
  bool has_P_{true};
};
