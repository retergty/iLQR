#pragma once
#include "Types.hpp"
#include "Dynamics.hpp"
#include "LinearApproximation.hpp"
enum class SensitivityIntegratorType
{
  EULER,
  RK2,
  RK4
};

template <typename Scalar, int XDimisions, int UDimisions>
class DynamicsDiscretizerBase
{
public:
  DynamicsDiscretizerBase() = default;
  virtual ~DynamicsDiscretizerBase() = default;

  /**
   * A function to computes the discrete approximation of the system's flowmap.
   * @param system : system to be discretized
   * @param t : starting time of the discretization interval
   * @param x : starting state x_{k}
   * @param u : input u_{k}, assumed constant over the entire interval
   * @param dt : interval duration
   * Returns x_{k+1}
   */
  virtual Vector<Scalar, XDimisions> discretize(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, const Scalar t, const Vector<Scalar, XDimisions> &x,
                                                const Vector<Scalar, UDimisions> &u, const Scalar dt) = 0;

  /**
   * A function to compute the linear approximation of the discretized system's flowmap.
   *
   * @param system : system to be discretized
   * @param t : starting time of the discretization interval
   * @param x : starting state x_{k}
   * @param u : input u_{k}, assumed constant over the entire interval
   * @param dt : interval duration
   * Returns an approximation of the form:
   *      x_{k+1} = A_{k} * dx_{k} + B_{k} * du_{k} + b_{k}
   */
  virtual VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>
  sensitivityDiscretize(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, const Scalar t, const Vector<Scalar, XDimisions> &x,
                        const Vector<Scalar, UDimisions> &u, const Scalar dt) = 0;
};

// Uses an Forward euler discretization.
template <typename Scalar, int XDimisions, int UDimisions>
class EulerDynamicsDiscretizer : public DynamicsDiscretizerBase<Scalar, XDimisions, UDimisions>
{
public:
  Vector<Scalar, XDimisions> discretize(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, const Scalar t, const Vector<Scalar, XDimisions> &x,
                                        const Vector<Scalar, UDimisions> &u, const Scalar dt) override
  {
    Vector<Scalar, XDimisions> tmp = system.computeFlowMap(t, x, u);
    tmp = x + dt * tmp;
    return tmp;
  }
  VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>
  sensitivityDiscretize(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, const Scalar t, const Vector<Scalar, XDimisions> &x,
                        const Vector<Scalar, UDimisions> &u, const Scalar dt)
  {
    // x_{k+1} = A_{k} * dx_{k} + B_{k} * du_{k} + b_{k}
    // A_{k} = Id + dt * dfdx
    // B_{k} = dt * dfdu
    // b_{k} = x_{n} + dt * f(x_{n},u_{n})
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> continuousApproximation = system.linearApproximation(t, x, u);
    continuousApproximation.dfdx *= dt;
    continuousApproximation.dfdx.diagonal().array() += 1.0; // plus Identity()
    continuousApproximation.dfdu *= dt;
    continuousApproximation.f = x + dt * continuousApproximation.f;
    return continuousApproximation;
  }
};

// Uses an Runge-Kutta 2nd order discretization
template <typename Scalar, int XDimisions, int UDimisions>
class EK2DynamicsDiscretizer : public DynamicsDiscretizerBase<Scalar, XDimisions, UDimisions>
{
public:
  Vector<Scalar, XDimisions> discretize(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, const Scalar t, const Vector<Scalar, XDimisions> &x,
                                        const Vector<Scalar, UDimisions> &u, const Scalar dt) override
  {
    const Scalar dt_halve = dt / 2.0;

    // System evaluations
    const Vector<Scalar, XDimisions> k1 = system.computeFlowMap(t, x, u);

    Vector<Scalar, XDimisions> tmp = x + dt * k1;
    const Vector<Scalar, XDimisions> k2 = system.computeFlowMap(t + dt, tmp, u);

    tmp = x + dt_halve * k1 + dt_halve * k2;
    return tmp;
  }
  VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>
  sensitivityDiscretize(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, const Scalar t, const Vector<Scalar, XDimisions> &x,
                        const Vector<Scalar, UDimisions> &u, const Scalar dt)
  {
    const Scalar dt_halve = dt / 2.0;

    // System evaluations
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> k1 = system.linearApproximation(t, x, u);
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> k2 = system.linearApproximation(t + dt, x + dt * k1.f, u);

    // Input sensitivity \dot{Su} = dfdx(t) Su + dfdu(t), with Su(0) = Zero()
    // Re-use memory from k.dfdu as dkduk
    // dk1duk = k1.dfdu
    k2.dfdu += dt * k2.dfdx * k1.dfdu;

    // State sensitivity \dot{Sx} = dfdx(t) Sx, with Sx(0) = Identity()
    // Re-use memory from k.dfdx as dkdxk
    // dk1dxk = k1.dfdx;
    k2.dfdx += dt * k2.dfdx * k1.dfdx; // need one temporary to avoid alias

    // Assemble discrete approximation
    // Re-use k1 to collect the result
    k1.dfdx = dt_halve * k1.dfdx + dt_halve * k2.dfdx;
    k1.dfdx.diagonal().array() += 1.0; // plus Identity()
    k1.dfdu = dt_halve * k1.dfdu + dt_halve * k2.dfdu;
    k1.f = x + dt_halve * k1.f + dt_halve * k2.f;
    return k1;
  }
};

// Uses an Runge-Kutta 4th order discretization.
template <typename Scalar, int XDimisions, int UDimisions>
class EK4DynamicsDiscretizer : public DynamicsDiscretizerBase<Scalar, XDimisions, UDimisions>
{
public:
  Vector<Scalar, XDimisions> discretize(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, const Scalar t, const Vector<Scalar, XDimisions> &x,
                                        const Vector<Scalar, UDimisions> &u, const Scalar dt) override
  {
    const Scalar dt_halve = dt / 2.0;
    const Scalar dt_sixth = dt / 6.0;
    const Scalar dt_third = dt / 3.0;

    // System evaluations
    const Vector<Scalar, XDimisions> k1 = system.computeFlowMap(t, x, u);
    Vector<Scalar, XDimisions> tmp = x + dt_halve * k1;
    const Vector<Scalar, XDimisions> k2 = system.computeFlowMap(t + dt_halve, tmp, u);
    tmp = x + dt_halve * k2;
    const Vector<Scalar, XDimisions> k3 = system.computeFlowMap(t + dt_halve, tmp, u);
    tmp = x + dt * k3;
    const Vector<Scalar, XDimisions> k4 = system.computeFlowMap(t + dt, tmp, u);

    tmp = x + dt_sixth * k1 + dt_third * k2 + dt_third * k3 + dt_sixth * k4;
    return tmp;
  }
  VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>
  sensitivityDiscretize(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, const Scalar t, const Vector<Scalar, XDimisions> &x,
                        const Vector<Scalar, UDimisions> &u, const Scalar dt)
  {
    const Scalar dt_halve = dt / 2.0;
    const Scalar dt_sixth = dt / 6.0;
    const Scalar dt_third = dt / 3.0;

    // System evaluations
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> k1 = system.linearApproximation(t, x, u);
    Vector<Scalar, XDimisions> tmpV = x + dt_halve * k1.f;
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> k2 = system.linearApproximation(t + dt_halve, tmpV, u);
    tmpV = x + dt_halve * k2.f;
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> k3 = system.linearApproximation(t + dt_halve, tmpV, u);
    tmpV = x + dt * k3.f;
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> k4 = system.linearApproximation(t + dt, tmpV, u);

    // Input sensitivity \dot{Su} = dfdx(t) Su + dfdu(t), with Su(0) = Zero()
    // Re-use memory from k.dfdu as dkduk
    // dk1duk = k1.dfdu
    k2.dfdu += dt_halve * k2.dfdx * k1.dfdu;
    k3.dfdu += dt_halve * k3.dfdx * k2.dfdu;
    k4.dfdu += dt * k4.dfdx * k3.dfdu;

    // State sensitivity \dot{Sx} = dfdx(t) Sx, with Sx(0) = Identity()
    // Re-use memory from k.dfdx as dkdxk
    // dk1dxk = k1.dfdx;
    Matrix<Scalar, XDimisions, XDimisions> tmp = dt_halve * k2.dfdx * k1.dfdx; // need one temporary to avoid alias
    k2.dfdx += tmp;
    tmp = dt_halve * k3.dfdx * k2.dfdx;
    k3.dfdx += tmp;
    tmp = dt * k4.dfdx * k3.dfdx;
    k4.dfdx += tmp;

    // Assemble discrete approximation
    // Re-use k1 to collect the result
    k1.dfdx = dt_sixth * k1.dfdx + dt_third * k2.dfdx + dt_third * k3.dfdx + dt_sixth * k4.dfdx;
    k1.dfdx += 1.0; // plus Identity()
    k1.dfdu = dt_sixth * k1.dfdu + dt_third * k2.dfdu + dt_third * k3.dfdu + dt_sixth * k4.dfdu;
    k1.f = x + dt_sixth * k1.f + dt_third * k2.f + dt_third * k3.f + dt_sixth * k4.f;
    return k1;
  }
};