/**
 * @file SensitivityIntegrator.hpp
 * @brief 动力学离散化与灵敏度：将连续时间动力学离散为
 * x_{k+1}=A*dx+B*du+b，支持欧拉/RK2/RK4。
 */
#pragma once
#include "Approximation/LinearApproximation.hpp"
#include "Types.hpp"

template <typename Scalar, int XDim, int UDim>
class SystemDynamicsBase;

/** @brief 灵敏度积分器类型：欧拉、RK2、RK4。 */
enum class SensitivityIntegratorType { EULER, RK2, RK4 };

/**
 * @brief 动力学离散化基类：由 (t,x,u,dt) 计算离散状态及线性近似 (A,B,b)。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int XDim, int UDim>
class DynamicsDiscretizerBase {
 public:
  DynamicsDiscretizerBase() = default;
  virtual ~DynamicsDiscretizerBase() = default;

  /**
   * @brief 计算离散化流映射 x_{k+1}。
   * @param [in] system 待离散的系统动力学。
   * @param [in] t 区间起始时间。
   * @param [in] x 起始状态 x_k。
   * @param [in] u 输入 u_k（区间内常值）。
   * @param [in] dt 区间长度。
   * @return 下一状态 x_{k+1}。
   */
  virtual Vector<Scalar, XDim> discretize(
      SystemDynamicsBase<Scalar, XDim, UDim>& system, const Scalar t,
      const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u,
      const Scalar dt) = 0;

  /**
   * @brief 计算离散化流映射的线性近似 x_{k+1} ≈ A*dx + B*du + b。
   * @param [in] system 待离散的系统动力学。
   * @param [in] t 区间起始时间。
   * @param [in] x 起始状态。
   * @param [in] u 输入。
   * @param [in] dt 区间长度。
   * @return 线性近似 (f, dfdx, dfdu)。
   */
  virtual VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
  sensitivityDiscretize(SystemDynamicsBase<Scalar, XDim, UDim>& system,
                        const Scalar t, const Vector<Scalar, XDim>& x,
                        const Vector<Scalar, UDim>& u, const Scalar dt) = 0;
};

/** @brief 前向欧拉离散化。 */
template <typename Scalar, int XDim, int UDim>
class EulerDynamicsDiscretizer
    : public DynamicsDiscretizerBase<Scalar, XDim, UDim> {
 public:
  Vector<Scalar, XDim> discretize(
      SystemDynamicsBase<Scalar, XDim, UDim>& system, const Scalar t,
      const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u,
      const Scalar dt) override {
    Vector<Scalar, XDim> tmp = system.computeFlowMap(t, x, u);
    tmp = x + dt * tmp;
    return tmp;
  }
  VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
  sensitivityDiscretize(SystemDynamicsBase<Scalar, XDim, UDim>& system,
                        const Scalar t, const Vector<Scalar, XDim>& x,
                        const Vector<Scalar, UDim>& u,
                        const Scalar dt) override {
    // x_{k+1} = A_{k} * dx_{k} + B_{k} * du_{k} + b_{k}
    // A_{k} = Id + dt * dfdx。
    // B_{k} = dt * dfdu。
    // b_{k} = x_{n} + dt * f(x_{n},u_{n})
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
        continuousApproximation = system.linearApproximation(t, x, u);
    continuousApproximation.dfdx *= dt;
    continuousApproximation.dfdx +=
        Matrix<Scalar, XDim, XDim>::Identity();  // 加上 Identity()。
    continuousApproximation.dfdu *= dt;
    continuousApproximation.f = x + dt * continuousApproximation.f;
    return continuousApproximation;
  }
};

// 使用二阶 Runge-Kutta 离散化。
template <typename Scalar, int XDim, int UDim>
class EK2DynamicsDiscretizer
    : public DynamicsDiscretizerBase<Scalar, XDim, UDim> {
 public:
  Vector<Scalar, XDim> discretize(
      SystemDynamicsBase<Scalar, XDim, UDim>& system, const Scalar t,
      const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u,
      const Scalar dt) override {
    const Scalar dt_halve = dt / 2.0;

    // 系统评估
    const Vector<Scalar, XDim> k1 = system.computeFlowMap(t, x, u);

    Vector<Scalar, XDim> tmp = x + dt * k1;
    const Vector<Scalar, XDim> k2 = system.computeFlowMap(t + dt, tmp, u);

    tmp = x + dt_halve * k1 + dt_halve * k2;
    return tmp;
  }
  VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
  sensitivityDiscretize(SystemDynamicsBase<Scalar, XDim, UDim>& system,
                        const Scalar t, const Vector<Scalar, XDim>& x,
                        const Vector<Scalar, UDim>& u,
                        const Scalar dt) override {
    const Scalar dt_halve = dt / 2.0;

    // 系统评估
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> k1 =
        system.linearApproximation(t, x, u);
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> k2 =
        system.linearApproximation(t + dt, x + dt * k1.f, u);

    // 输入灵敏度 \dot{Su} = dfdx(t) Su + dfdu(t)，且 Su(0) = Zero()。
    // 复用 k.dfdu 的内存作为 dkduk。
    // dk1duk = k1.dfdu
    k2.dfdu += dt * k2.dfdx * k1.dfdu;

    // 状态灵敏度 \dot{Sx} = dfdx(t) Sx，且 Sx(0) = Identity()。
    // 复用 k.dfdx 的内存作为 dkdxk。
    // dk1dxk = k1.dfdx;
    k2.dfdx += dt * k2.dfdx * k1.dfdx;  // 需要一个临时变量以避免别名。

    // 组装离散近似。
    // 复用 k1 收集结果。
    k1.dfdx = dt_halve * k1.dfdx + dt_halve * k2.dfdx;
    k1.dfdx += Matrix<Scalar, XDim, XDim>::Identity();  // 加上 Identity()。
    k1.dfdu = dt_halve * k1.dfdu + dt_halve * k2.dfdu;
    k1.f = x + dt_halve * k1.f + dt_halve * k2.f;
    return k1;
  }
};

// 使用四阶 Runge-Kutta 离散化。
template <typename Scalar, int XDim, int UDim>
class EK4DynamicsDiscretizer
    : public DynamicsDiscretizerBase<Scalar, XDim, UDim> {
 public:
  Vector<Scalar, XDim> discretize(
      SystemDynamicsBase<Scalar, XDim, UDim>& system, const Scalar t,
      const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u,
      const Scalar dt) override {
    const Scalar dt_halve = dt / 2.0;
    const Scalar dt_sixth = dt / 6.0;
    const Scalar dt_third = dt / 3.0;

    // 系统评估
    const Vector<Scalar, XDim> k1 = system.computeFlowMap(t, x, u);
    Vector<Scalar, XDim> tmp = x + dt_halve * k1;
    const Vector<Scalar, XDim> k2 = system.computeFlowMap(t + dt_halve, tmp, u);
    tmp = x + dt_halve * k2;
    const Vector<Scalar, XDim> k3 = system.computeFlowMap(t + dt_halve, tmp, u);
    tmp = x + dt * k3;
    const Vector<Scalar, XDim> k4 = system.computeFlowMap(t + dt, tmp, u);

    tmp = x + dt_sixth * k1 + dt_third * k2 + dt_third * k3 + dt_sixth * k4;
    return tmp;
  }
  VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
  sensitivityDiscretize(SystemDynamicsBase<Scalar, XDim, UDim>& system,
                        const Scalar t, const Vector<Scalar, XDim>& x,
                        const Vector<Scalar, UDim>& u,
                        const Scalar dt) override {
    const Scalar dt_halve = dt / 2.0;
    const Scalar dt_sixth = dt / 6.0;
    const Scalar dt_third = dt / 3.0;

    // 系统评估
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> k1 =
        system.linearApproximation(t, x, u);
    Vector<Scalar, XDim> tmpV = x + dt_halve * k1.f;
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> k2 =
        system.linearApproximation(t + dt_halve, tmpV, u);
    tmpV = x + dt_halve * k2.f;
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> k3 =
        system.linearApproximation(t + dt_halve, tmpV, u);
    tmpV = x + dt * k3.f;
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> k4 =
        system.linearApproximation(t + dt, tmpV, u);

    // 输入灵敏度 \dot{Su} = dfdx(t) Su + dfdu(t)，且 Su(0) = Zero()。
    // 复用 k.dfdu 的内存作为 dkduk。
    // dk1duk = k1.dfdu
    k2.dfdu += dt_halve * k2.dfdx * k1.dfdu;
    k3.dfdu += dt_halve * k3.dfdx * k2.dfdu;
    k4.dfdu += dt * k4.dfdx * k3.dfdu;

    // 状态灵敏度 \dot{Sx} = dfdx(t) Sx，且 Sx(0) = Identity()。
    // 复用 k.dfdx 的内存作为 dkdxk。
    // dk1dxk = k1.dfdx;
    Matrix<Scalar, XDim, XDim> tmp =
        dt_halve * k2.dfdx * k1.dfdx;  // 需要一个临时变量以避免别名。
    k2.dfdx += tmp;
    tmp = dt_halve * k3.dfdx * k2.dfdx;
    k3.dfdx += tmp;
    tmp = dt * k4.dfdx * k3.dfdx;
    k4.dfdx += tmp;

    // 组装离散近似。
    // 复用 k1 收集结果。
    k1.dfdx = dt_sixth * k1.dfdx + dt_third * k2.dfdx + dt_third * k3.dfdx +
              dt_sixth * k4.dfdx;
    k1.dfdx += Matrix<Scalar, XDim, XDim>::Identity();  // 加上 Identity()。
    k1.dfdu = dt_sixth * k1.dfdu + dt_third * k2.dfdu + dt_third * k3.dfdu +
              dt_sixth * k4.dfdu;
    k1.f = x + dt_sixth * k1.f + dt_third * k2.f + dt_third * k3.f +
           dt_sixth * k4.f;
    return k1;
  }
};