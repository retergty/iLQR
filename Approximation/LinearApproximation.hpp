/**
 * @file LinearApproximation.hpp
 * @brief 线性近似：标量/向量函数的一阶近似结构（f(x,u) ≈ dfdx' dx + dfdu' du +
 * f）。
 */
#pragma once
#include "Matrix/Types.hpp"

/**
 * @brief 标量函数线性近似：f(x,u) = dfdx' dx + dfdu' du + f。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int XDim, int UDim>
struct ScalarFunctionLinearApproximation {
  /** @brief 对状态的一阶导数。 */
  Vector<Scalar, XDim> dfdx;
  /** @brief 对输入的一阶导数。 */
  Vector<Scalar, UDim> dfdu;
  /** @brief 常数项。 */
  Scalar f{0};

  /** @brief 默认构造。 */
  ScalarFunctionLinearApproximation() = default;

  /** @brief 复合加法赋值。 */
  ScalarFunctionLinearApproximation& operator+=(
      const ScalarFunctionLinearApproximation& rhs) {
    dfdx += rhs.dfdx;
    dfdu += rhs.dfdu;
    f += rhs.f;
    return *this;
  }

  /** @brief 复合标量乘法赋值。 */
  ScalarFunctionLinearApproximation& operator*=(Scalar s) {
    dfdx *= s;
    dfdu *= s;
    f *= s;
    return *this;
  }

  /** @brief 将各系数置零。 */
  ScalarFunctionLinearApproximation& setZero() {
    dfdx.setZero();
    dfdu.setZero();
    f = 0;
    return *this;
  }

  /** @brief 返回零初始化的近似对象。 */
  static ScalarFunctionLinearApproximation Zero() {
    ScalarFunctionLinearApproximation f;
    f.setZero();
    return f;
  }
};

/** @brief 无输入时的标量函数线性近似特化（仅 dfdx 与 f）。 */
template <typename Scalar, int XDim>
struct ScalarFunctionLinearApproximation<Scalar, XDim, 0> {
  /** @brief 对状态的一阶导数。 */
  Vector<Scalar, XDim> dfdx;
  /** @brief 常数项。 */
  Scalar f{0};

  /** @brief 默认构造。 */
  ScalarFunctionLinearApproximation() = default;

  /** @brief 复合加法赋值。 */
  ScalarFunctionLinearApproximation& operator+=(
      const ScalarFunctionLinearApproximation& rhs) {
    dfdx += rhs.dfdx;
    f += rhs.f;
    return *this;
  }

  /** @brief 复合标量乘法赋值。 */
  ScalarFunctionLinearApproximation& operator*=(Scalar s) {
    dfdx *= s;
    f *= s;
    return *this;
  }

  /** @brief 将各系数置零。 */
  ScalarFunctionLinearApproximation& setZero() {
    dfdx.setZero();
    f = 0;
    return *this;
  }

  /** @brief 返回零初始化的近似对象。 */
  static ScalarFunctionLinearApproximation Zero() {
    ScalarFunctionLinearApproximation f;
    f.setZero();
    return f;
  }
};

/**
 * @brief 向量函数线性近似：f(x,u) = dfdx*dx + dfdu*du + f。
 * @tparam Scalar 标量类型。
 * @tparam FDimisions 函数输出维度。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int FDimisions, int XDim, int UDim>
struct VectorFunctionLinearApproximation {
  /** @brief 对状态的雅可比。 */
  Matrix<Scalar, FDimisions, XDim> dfdx;
  /** @brief 对输入的雅可比。 */
  Matrix<Scalar, FDimisions, UDim> dfdu;
  /** @brief 常数项向量。 */
  Vector<Scalar, FDimisions> f;

  /** @brief 默认构造。 */
  VectorFunctionLinearApproximation() = default;

  /** @brief 将各系数置零。 */
  VectorFunctionLinearApproximation& setZero() {
    dfdx.setZero();
    dfdu.setZero();
    f.setZero();
    return *this;
  }

  /** @brief 返回零初始化的近似对象。 */
  static VectorFunctionLinearApproximation Zero() {
    VectorFunctionLinearApproximation f;
    f.setZero();
    return f;
  }
};

/** @brief 无输入时的向量函数线性近似特化。 */
template <typename Scalar, int FDimisions, int XDim>
struct VectorFunctionLinearApproximation<Scalar, FDimisions, XDim, 0> {
  /** @brief 对状态的雅可比。 */
  Matrix<Scalar, FDimisions, XDim> dfdx;
  /** @brief 常数项向量。 */
  Vector<Scalar, FDimisions> f;

  /** @brief 默认构造。 */
  VectorFunctionLinearApproximation() = default;

  /** @brief 将各系数置零。 */
  VectorFunctionLinearApproximation& setZero() {
    dfdx.setZero();
    f.setZero();
    return *this;
  }

  /** @brief 返回零初始化的近似对象。 */
  static VectorFunctionLinearApproximation Zero() {
    VectorFunctionLinearApproximation f;
    f.setZero();
    return f;
  }
};

template <typename Scalar, int XDim, int UDim>
ScalarFunctionLinearApproximation<Scalar, XDim, UDim>& operator+=(
    ScalarFunctionLinearApproximation<Scalar, XDim, UDim>& lhs,
    ScalarFunctionLinearApproximation<Scalar, XDim, 0>& rhs) {
  lhs.f += rhs.f;
  lhs.dfdx += rhs.dfdx;
  return lhs;
}