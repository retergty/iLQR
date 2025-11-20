#pragma once
#include <array>
#include "Types.hpp"

/**
 * Defines the linear approximation of a scalar function
 * f(x,u) = dfdx' dx + dfdu' du + f
 */
template <typename Scalar, int XDimisions, int UDimisions>
struct ScalarFunctionLinearApproximation
{
  /** First derivative w.r.t state */
  Vector<Scalar, XDimisions> dfdx;
  /** First derivative w.r.t input */
  Vector<Scalar, UDimisions> dfdu;
  /** Constant term */
  Scalar f{ 0 };

  /** Default constructor */
  ScalarFunctionLinearApproximation() = default;

  /** Compound addition assignment operator */
  ScalarFunctionLinearApproximation& operator+=(const ScalarFunctionLinearApproximation& rhs)
  {
    dfdx += rhs.dfdx;
    dfdu += rhs.dfdu;
    f += rhs.f;
    return *this;
  }

  /** Compound scalar multiplication and assignment operator */
  ScalarFunctionLinearApproximation& operator*=(Scalar s)
  {
    dfdx *= s;
    dfdu *= s;
    f *= s;
    return *this;
  }

  /**
   * Resizes the members to the given size, and sets all coefficients to zero.
   * @param[in] nx State dimension
   * @param[in] nu Input dimension (Pass nu = -1 for no inputs)
   */
  ScalarFunctionLinearApproximation& setZero()
  {
    dfdx.setZero();
    dfdu.setZero();
    f = 0;
    return *this;
  }

  /**
   * Factory function with zero initialization
   * @param[in] nx State dimension
   * @param[in] nu Input dimension (Pass nu = -1 for no inputs)
   * @return Zero initialized object of given size.
   */
  static ScalarFunctionLinearApproximation Zero()
  {
    ScalarFunctionLinearApproximation f;
    f.setZero();
    return f;
  }
};

// no input
template <typename Scalar, int XDimisions>
struct ScalarFunctionLinearApproximation<Scalar, XDimisions, 0>
{
  /** First derivative w.r.t state */
  Vector<Scalar, XDimisions> dfdx;
  /** Constant term */
  Scalar f{ 0 };

  /** Default constructor */
  ScalarFunctionLinearApproximation() = default;

  /** Compound addition assignment operator */
  ScalarFunctionLinearApproximation& operator+=(const ScalarFunctionLinearApproximation& rhs)
  {
    dfdx += rhs.dfdx;
    f += rhs.f;
    return *this;
  }

  /** Compound scalar multiplication and assignment operator */
  ScalarFunctionLinearApproximation& operator*=(Scalar s)
  {
    dfdx *= s;
    f *= s;
    return *this;
  }

  /**
   * sets all coefficients to zero.
   */
  ScalarFunctionLinearApproximation& setZero()
  {
    dfdx.setZero();
    f = 0;
  }

  /**
   * Factory function with zero initialization
   * @return Zero initialized object.
   */
  static ScalarFunctionLinearApproximation Zero()
  {
    ScalarFunctionLinearApproximation f;
    f.setZero();
    return f;
  }
};

// no input
template <typename Scalar, int FDimisions, int XDimisions, int UDimisions>
struct VectorFunctionLinearApproximation
{
  /** Derivative w.r.t state */
  Matrix<Scalar, FDimisions, XDimisions> dfdx;
  /** Derivative w.r.t input */
  Matrix<Scalar, FDimisions, UDimisions> dfdu;
  /** Constant term */
  Vector<Scalar, FDimisions> f;

  /** Default constructor */
  VectorFunctionLinearApproximation() = default;

  /**
   * sets all coefficients to zero.
   */
  VectorFunctionLinearApproximation& setZero()
  {
    dfdx.setZero();
    dfdu.setZero();
    f.setZero();
    return *this;
  }

  /**
   * Factory function with zero initialization
   * @return Zero initialized object.
   */
  static VectorFunctionLinearApproximation Zero()
  {
    VectorFunctionLinearApproximation f;
    f.setZero();
    return f;
  }
};

/**
 * Defines the linear model of a vector-valued function
 * f(x,u) = \nabla_x^T dx + \nabla_u^T du + f
 */
template <typename Scalar, int FDimisions, int XDimisions>
struct VectorFunctionLinearApproximation<Scalar, FDimisions, XDimisions, 0>
{
  /** Derivative w.r.t state */
  Matrix<Scalar, FDimisions, XDimisions> dfdx;
  /** Constant term */
  Vector<Scalar, FDimisions> f;

  /** Default constructor */
  VectorFunctionLinearApproximation() = default;

  /**
   * sets all coefficients to zero.
   */
  VectorFunctionLinearApproximation& setZero()
  {
    dfdx.setZero();
    f.setZero();
    return *this;
  }

  /**
   * Factory function with zero initialization
   * @return Zero initialized object.
   */
  static VectorFunctionLinearApproximation Zero()
  {
    VectorFunctionLinearApproximation f;
    f.setZero();
    return f;
  }
};

template<typename Scalar, int XDimisions, int UDimisions>
ScalarFunctionLinearApproximation<Scalar, XDimisions, UDimisions>& operator+=(ScalarFunctionLinearApproximation<Scalar, XDimisions, UDimisions>& lhs, ScalarFunctionLinearApproximation<Scalar, XDimisions, 0>& rhs)
{
  lhs.f += rhs.f;
  lhs.dfdx += rhs.dfdx;
  return lhs;
}