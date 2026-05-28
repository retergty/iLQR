/**
 * @file DiscreteLinearSystemDynamics.hpp
 * @brief 线性离散时不变系统动力学：
 * \f$ x_{k+1} = A x_k + B u_k + b \f$。
 */
#pragma once

#include "Dynamics/DiscreteSystemDynamicsBase.hpp"

/**
 * @brief 线性离散时不变系统动力学。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class DiscreteLinearSystemDynamics
    : public DiscreteSystemDynamicsBase<Scalar, XDim, UDim> {
 public:
  using LinearApproximation_t =
      typename DiscreteSystemDynamicsBase<Scalar, XDim,
                                          UDim>::LinearApproximation_t;

  /**
   * @brief 用状态矩阵 A 与控制矩阵 B 构造无偏置离散系统。
   * @param [in] A 离散状态转移矩阵。
   * @param [in] B 离散控制矩阵。
   */
  DiscreteLinearSystemDynamics(const Matrix<Scalar, XDim, XDim>& A,
                               const Matrix<Scalar, XDim, UDim>& B)
      : A_(A), B_(B) {
    bias_.setZero();
  }

  /**
   * @brief 用状态矩阵 A、控制矩阵 B 与仿射偏置 b 构造离散系统。
   * @param [in] A 离散状态转移矩阵。
   * @param [in] B 离散控制矩阵。
   * @param [in] bias 仿射偏置。
   */
  DiscreteLinearSystemDynamics(const Matrix<Scalar, XDim, XDim>& A,
                               const Matrix<Scalar, XDim, UDim>& B,
                               const Vector<Scalar, XDim>& bias)
      : A_(A), B_(B), bias_(bias) {}

  /** @brief 析构函数。 */
  ~DiscreteLinearSystemDynamics() override = default;

  /** @brief 计算离散状态转移：\f$ x_{k+1}=A x_k+B u_k+b \f$。 */
  Vector<Scalar, XDim> computeMap(Scalar t, const Vector<Scalar, XDim>& x,
                                  const Vector<Scalar, UDim>& u,
                                  Scalar dt) override {
    (void)t;
    (void)dt;
    Vector<Scalar, XDim> nextState = A_ * x;
    nextState += B_ * u;
    nextState += bias_;
    return nextState;
  }

  /** @brief 返回完整离散线性近似：dfdx=A, dfdu=B, f=A*x+B*u+b。 */
  LinearApproximation_t linearApproximation(Scalar t,
                                            const Vector<Scalar, XDim>& x,
                                            const Vector<Scalar, UDim>& u,
                                            Scalar dt) override {
    LinearApproximation_t approximation;
    approximation.f = computeMap(t, x, u, dt);
    approximation.dfdx = A_;
    approximation.dfdu = B_;
    return approximation;
  }

  /** @brief 返回离散偏差动力学近似：dfdx=A, dfdu=B, f=0。 */
  LinearApproximation_t deviationLinearApproximation(
      Scalar t, const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u,
      Scalar dt) override {
    (void)t;
    (void)x;
    (void)u;
    (void)dt;
    LinearApproximation_t approximation;
    approximation.f.setZero();
    approximation.dfdx = A_;
    approximation.dfdu = B_;
    return approximation;
  }

 protected:
  /** @brief 拷贝构造（供子类使用）。 */
  DiscreteLinearSystemDynamics(const DiscreteLinearSystemDynamics& other) =
      default;

  /** @brief 离散状态转移矩阵。 */
  Matrix<Scalar, XDim, XDim> A_;
  /** @brief 离散控制矩阵。 */
  Matrix<Scalar, XDim, UDim> B_;
  /** @brief 离散仿射偏置。 */
  Vector<Scalar, XDim> bias_;
};
