#pragma once
#include <array>
#include "Types.hpp"

/**
 * Defines the quadratic approximation of a scalar function
 * f(x,u) = 1/2 dx' dfdxx dx + du' dfdux dx + 1/2 du' dfduu du + dfdx' dx + dfdu' du + f
 */
template <typename Scalar, int XDimisions, int UDimisions>
struct ScalarFunctionQuadraticApproximation
{
  /** Second derivative w.r.t state */
  Matrix<Scalar, XDimisions, XDimisions> dfdxx;
  /** Second derivative w.r.t input (lhs) and state (rhs) */
  Matrix<Scalar, UDimisions, XDimisions> dfdux;
  /** Second derivative w.r.t input */
  Matrix<Scalar, UDimisions, UDimisions> dfduu;
  /** First derivative w.r.t state */
  Vector<Scalar, XDimisions> dfdx;
  /** First derivative w.r.t input */
  Vector<Scalar, UDimisions> dfdu;
  /** Constant term */
  Scalar f{ 0 };

  /** Default constructor */
  ScalarFunctionQuadraticApproximation() = default;

  ScalarFunctionQuadraticApproximation(const ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>& rhs) : dfdxx(rhs.dfdxx), dfdx(rhs.dfdx), f(rhs.f)
  {
    dfdux.setZero();
    dfduu.setZero();
    dfdu.setZero();
  }

  /** Compound addition assignment operator */
  ScalarFunctionQuadraticApproximation& operator+=(const ScalarFunctionQuadraticApproximation& rhs)
  {
    dfdxx += rhs.dfdxx;
    dfdux += rhs.dfdux;
    dfduu += rhs.dfduu;
    dfdx += rhs.dfdx;
    dfdu += rhs.dfdu;
    f += rhs.f;

    return *this;
  }

  /** Compound scalar multiplication and assignment operator */
  ScalarFunctionQuadraticApproximation& operator*=(Scalar s)
  {
    dfdxx *= s;
    dfdux *= s;
    dfduu *= s;
    dfdx *= s;
    dfdu *= s;
    f *= s;

    return *this;
  }

  /**
   * sets all coefficients to zero.
   */
  ScalarFunctionQuadraticApproximation& setZero()
  {
    dfdxx.setZero();
    dfduu.setZero();
    dfdux.setZero();
    dfdu.setZero();
    dfdx.setZero();
    f = 0;
    return *this;
  }

  /**
   * Factory function with zero initialization
   * @return Zero initialized object.
   */
  static ScalarFunctionQuadraticApproximation Zero()
  {
    ScalarFunctionQuadraticApproximation f;
    f.setZero();

    return f;
  }
};

/** no input
 */
template <typename Scalar, int XDimisions>
struct ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
{
  /** Second derivative w.r.t state */
  Matrix<Scalar, XDimisions, XDimisions> dfdxx;
  /** First derivative w.r.t state */
  Vector<Scalar, XDimisions> dfdx;
  /** Constant term */
  Scalar f{ 0 };

  /** Default constructor */
  ScalarFunctionQuadraticApproximation() = default;

  /** Compound addition assignment operator */
  ScalarFunctionQuadraticApproximation& operator+=(const ScalarFunctionQuadraticApproximation& rhs)
  {
    dfdxx += rhs.dfdxx;
    dfdx += rhs.dfdx;
    f += rhs.f;

    return *this;
  }

  /** Compound scalar multiplication and assignment operator */
  ScalarFunctionQuadraticApproximation& operator*=(Scalar s)
  {
    dfdxx *= s;
    dfdx *= s;
    f *= s;

    return *this;
  }

  /**
   * Resizes the members to the given size, and sets all coefficients to zero.
   */
  ScalarFunctionQuadraticApproximation& setZero()
  {
    dfdxx.setZero();
    dfdx.setZero();
    f = 0;
    return *this;
  }

  /**
   * Factory function with zero initialization
   * @return Zero initialized object.
   */
  static ScalarFunctionQuadraticApproximation Zero()
  {
    ScalarFunctionQuadraticApproximation f;
    f.setZero();

    return f;
  }
};

template <typename Scalar, int XDimisions, int UDimisions>
inline ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> operator*(ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> lhs, Scalar scalar) {
  return lhs *= scalar;
}
template <typename Scalar, int XDimisions, int UDimisions>
inline ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> operator*(Scalar scalar, ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> rhs) {
  return rhs *= scalar;
}

template<typename Scalar, int XDimisions, int UDimisions>
ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>& operator+=(ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>& lhs, ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>& rhs)
{
  lhs.f += rhs.f;
  lhs.dfdx += rhs.dfdx;
  lhs.dfdxx += rhs.dfdxx;
  return lhs;
}