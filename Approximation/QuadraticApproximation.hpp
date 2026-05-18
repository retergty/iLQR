/**
 * @file QuadraticApproximation.hpp
 * @brief 二次近似：标量函数的二阶近似结构（Hessian、一阶项与常数项）。
 */
#pragma once
#include "Types.hpp"
#include <array>

/**
 * @brief 标量函数二次近似：f = 1/2 dx' dfdxx dx + du' dfdux dx + 1/2 du' dfduu
 * du + dfdx' dx + dfdu' du + f。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int XDim, int UDim>
struct ScalarFunctionQuadraticApproximation {
  /** @brief 对状态的二阶导数（Hessian）。 */
  Matrix<Scalar, XDim, XDim> dfdxx;
  /** @brief 对输入-状态的混合二阶导数。 */
  Matrix<Scalar, UDim, XDim> dfdux;
  /** @brief 对输入的二阶导数。 */
  Matrix<Scalar, UDim, UDim> dfduu;
  /** @brief 对状态的一阶导数。 */
  Vector<Scalar, XDim> dfdx;
  /** @brief 对输入的一阶导数。 */
  Vector<Scalar, UDim> dfdu;
  /** @brief 常数项。 */
  Scalar f{0};

  /** @brief 默认构造。 */
  ScalarFunctionQuadraticApproximation() = default;

  ScalarFunctionQuadraticApproximation(
      const ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> &rhs)
      : dfdxx(rhs.dfdxx), dfdx(rhs.dfdx), f(rhs.f) {
    dfdux.setZero();
    dfduu.setZero();
    dfdu.setZero();
  }

  /** @brief 复合加法赋值。 */
  ScalarFunctionQuadraticApproximation &
  operator+=(const ScalarFunctionQuadraticApproximation &rhs) {
    dfdxx += rhs.dfdxx;
    dfdux += rhs.dfdux;
    dfduu += rhs.dfduu;
    dfdx += rhs.dfdx;
    dfdu += rhs.dfdu;
    f += rhs.f;

    return *this;
  }

  /** @brief 复合标量乘法赋值。 */
  ScalarFunctionQuadraticApproximation &operator*=(Scalar s) {
    dfdxx *= s;
    dfdux *= s;
    dfduu *= s;
    dfdx *= s;
    dfdu *= s;
    f *= s;

    return *this;
  }

  /** @brief 将各系数置零。 */
  ScalarFunctionQuadraticApproximation &setZero() {
    dfdxx.setZero();
    dfduu.setZero();
    dfdux.setZero();
    dfdu.setZero();
    dfdx.setZero();
    f = 0;
    return *this;
  }

  /** @brief 返回零初始化的近似对象。 */
  static ScalarFunctionQuadraticApproximation Zero() {
    ScalarFunctionQuadraticApproximation f;
    f.setZero();

    return f;
  }
};

/** @brief 无输入时的标量函数二次近似特化（仅 dfdxx、dfdx、f）。 */
template <typename Scalar, int XDim>
struct ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> {
  /** @brief 对状态的二阶导数。 */
  Matrix<Scalar, XDim, XDim> dfdxx;
  /** @brief 对状态的一阶导数。 */
  Vector<Scalar, XDim> dfdx;
  /** @brief 常数项。 */
  Scalar f{0};

  /** @brief 默认构造。 */
  ScalarFunctionQuadraticApproximation() = default;

  /** @brief 复合加法赋值。 */
  ScalarFunctionQuadraticApproximation &
  operator+=(const ScalarFunctionQuadraticApproximation &rhs) {
    dfdxx += rhs.dfdxx;
    dfdx += rhs.dfdx;
    f += rhs.f;

    return *this;
  }

  /** @brief 复合标量乘法赋值。 */
  ScalarFunctionQuadraticApproximation &operator*=(Scalar s) {
    dfdxx *= s;
    dfdx *= s;
    f *= s;

    return *this;
  }

  /** @brief 将各系数置零。 */
  ScalarFunctionQuadraticApproximation &setZero() {
    dfdxx.setZero();
    dfdx.setZero();
    f = 0;
    return *this;
  }

  /** @brief 返回零初始化的近似对象。 */
  static ScalarFunctionQuadraticApproximation Zero() {
    ScalarFunctionQuadraticApproximation f;
    f.setZero();

    return f;
  }
};

template <typename Scalar, int XDim, int UDim>
inline ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
operator*(ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> lhs,
          Scalar scalar) {
  return lhs *= scalar;
}
template <typename Scalar, int XDim, int UDim>
inline ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
operator*(Scalar scalar,
          ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> rhs) {
  return rhs *= scalar;
}

template <typename Scalar, int XDim, int UDim>
ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> &
operator+=(ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> &lhs,
           ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> &rhs) {
  lhs.f += rhs.f;
  lhs.dfdx += rhs.dfdx;
  lhs.dfdxx += rhs.dfdxx;
  return lhs;
}