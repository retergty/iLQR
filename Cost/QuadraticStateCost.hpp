/**
 * @file QuadraticStateCost.hpp
 * @brief 二次仅状态代价：l = 0.5 (x-x_ref)' Q (x-x_ref)，支持参考轨迹插值。
 */
#pragma once

#include "Cost/Cost.hpp"
#include "Misc/LinearInterpolation.hpp"

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
   * @param [in] cost_number 代价项唯一标识。
   */
  QuadraticStateCost(const Matrix<Scalar, XDim, XDim>& Q, int cost_number)
      : StateCost<Scalar, XDim, ArrayLength>(cost_number), Q_(Q) {};
  ~QuadraticStateCost() override = default;

  /** @brief 获取代价值 0.5 * (x-x_ref)' Q (x-x_ref)。 */
  Scalar getValue(Scalar time, const Vector<Scalar, XDim>& state,
                  const std::array<Scalar, ArrayLength>& timeTrajectory,
                  const std::array<Vector<Scalar, XDim>, ArrayLength>&
                      stateTrajectoy) final {
    const Vector<Scalar, XDim> xDeviation =
        getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    return Scalar(0.5) * xDeviation.dot(Q_ * xDeviation);
  }

  /** @brief 计算并累加代价的二次近似（dfdxx=Q, dfdx=Q*(x-x_ref),
   * f=0.5*(x-x_ref)'*dfdx）。 */
  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>& addAppro) final {
    const Vector<Scalar, XDim> xDeviation =
        getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    const Vector<Scalar, XDim> weightedStateDeviation = Q_ * xDeviation;
    addAppro.dfdxx += Q_;
    addAppro.dfdx += weightedStateDeviation;
    addAppro.f += Scalar(0.5) * xDeviation.dot(weightedStateDeviation);
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
   * @param [in] cost_number 代价项唯一标识。
   */
  QuadraticStateInputCost(const Matrix<Scalar, XDim, XDim>& Q,
                          const Matrix<Scalar, UDim, UDim>& R,
                          const Matrix<Scalar, UDim, XDim>& P, int cost_number)
      : StateInputCost<Scalar, XDim, UDim, ArrayLength>(cost_number),
        Q_(Q),
        R_(R),
        P_(P) {
    has_P_ = true;
  };

  /** 用 Q、R 和代价项编号构造无交叉项的二次状态-输入代价。 */
  QuadraticStateInputCost(const Matrix<Scalar, XDim, XDim>& Q,
                          const Matrix<Scalar, UDim, UDim>& R, int cost_number)
      : StateInputCost<Scalar, XDim, UDim, ArrayLength>(cost_number),
        Q_(Q),
        R_(R) {
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
      final {
    Vector<Scalar, XDim> stateDeviation =
        getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    Vector<Scalar, UDim> inputDeviation =
        getInputDeviation(time, input, timeTrajectory, inputTrajectory);

    if (has_P_) {
      return Scalar(0.5) * stateDeviation.dot(Q_ * stateDeviation) +
             Scalar(0.5) * inputDeviation.dot(R_ * inputDeviation) +
             inputDeviation.dot(P_ * stateDeviation);
    } else {
      return Scalar(0.5) * stateDeviation.dot(Q_ * stateDeviation) +
             Scalar(0.5) * inputDeviation.dot(R_ * inputDeviation);
    }
  }

  /** 计算并累加代价项二次近似。 */
  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory,
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& addAppro)
      final {
    (void)timeTrajectory;
    Vector<Scalar, XDim> stateDeviation =
        getStateDeviation(time, state, timeTrajectory, stateTrajectoy);
    Vector<Scalar, UDim> inputDeviation =
        getInputDeviation(time, input, timeTrajectory, inputTrajectory);
    const Vector<Scalar, XDim> weightedStateDeviation = Q_ * stateDeviation;
    const Vector<Scalar, UDim> weightedInputDeviation = R_ * inputDeviation;
    addAppro.dfdxx += Q_;
    addAppro.dfduu += R_;
    addAppro.dfdx += weightedStateDeviation;
    addAppro.dfdu += weightedInputDeviation;
    addAppro.f += Scalar(0.5) * stateDeviation.dot(weightedStateDeviation) +
                  Scalar(0.5) * inputDeviation.dot(weightedInputDeviation);

    if (has_P_) {
      const Vector<Scalar, UDim> pDeviation = P_ * stateDeviation;
      addAppro.f += inputDeviation.dot(pDeviation);
      addAppro.dfdu += pDeviation;
      addAppro.dfdx += P_.transpose() * inputDeviation;
      addAppro.dfdux += P_;
    }
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
