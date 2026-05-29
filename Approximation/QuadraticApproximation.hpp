/**
 * @file QuadraticApproximation.hpp
 * @brief
 * 二次近似：标量/向量函数的二阶近似结构（Hessian、Jacobian/梯度与常数项）。
 */
#pragma once
#include <array>

#include "matrix/Types.hpp"

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
      const ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>& rhs)
      : dfdxx(rhs.dfdxx), dfdx(rhs.dfdx), f(rhs.f) {
    dfdux.setZero();
    dfduu.setZero();
    dfdu.setZero();
  }

  /** @brief 复合加法赋值。 */
  ScalarFunctionQuadraticApproximation& operator+=(
      const ScalarFunctionQuadraticApproximation& rhs) {
    dfdxx += rhs.dfdxx;
    dfdux += rhs.dfdux;
    dfduu += rhs.dfduu;
    dfdx += rhs.dfdx;
    dfdu += rhs.dfdu;
    f += rhs.f;

    return *this;
  }

  /** @brief 复合标量乘法赋值。 */
  ScalarFunctionQuadraticApproximation& operator*=(Scalar s) {
    dfdxx *= s;
    dfdux *= s;
    dfduu *= s;
    dfdx *= s;
    dfdu *= s;
    f *= s;

    return *this;
  }

  /** @brief 将各系数置零。 */
  ScalarFunctionQuadraticApproximation& setZero() {
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
  ScalarFunctionQuadraticApproximation& operator+=(
      const ScalarFunctionQuadraticApproximation& rhs) {
    dfdxx += rhs.dfdxx;
    dfdx += rhs.dfdx;
    f += rhs.f;

    return *this;
  }

  /** @brief 复合标量乘法赋值。 */
  ScalarFunctionQuadraticApproximation& operator*=(Scalar s) {
    dfdxx *= s;
    dfdx *= s;
    f *= s;

    return *this;
  }

  /** @brief 将各系数置零。 */
  ScalarFunctionQuadraticApproximation& setZero() {
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
inline ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> operator*(
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> lhs,
    Scalar scalar) {
  return lhs *= scalar;
}
template <typename Scalar, int XDim, int UDim>
inline ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> operator*(
    Scalar scalar,
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> rhs) {
  return rhs *= scalar;
}

template <typename Scalar, int XDim, int UDim>
ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& operator+=(
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& lhs,
    const ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>& rhs) {
  lhs.f += rhs.f;
  lhs.dfdx += rhs.dfdx;
  lhs.dfdxx += rhs.dfdxx;
  return lhs;
}

/**
 * @brief 向量函数的二次近似。
 *
 * 对每个输出分量 i：
 * f_i(x,u) = 1/2 dx' dfdxx[i] dx
 *          + du' dfdux[i] dx
 *          + 1/2 du' dfduu[i] du
 *          + dfdx[i,:] dx
 *          + dfdu[i,:] du
 *          + f[i]
 */
template <typename Scalar, int FDim, int XDim, int UDim>
struct VectorFunctionQuadraticApproximation {
  /** @brief 每个输出分量对状态的 Hessian。 */
  std::array<Matrix<Scalar, XDim, XDim>, FDim> dfdxx;
  /** @brief 每个输出分量对输入-状态的混合二阶导数，尺寸为 UDim x XDim。 */
  std::array<Matrix<Scalar, UDim, XDim>, FDim> dfdux;
  /** @brief 每个输出分量对输入的 Hessian。 */
  std::array<Matrix<Scalar, UDim, UDim>, FDim> dfduu;
  /** @brief 对状态的 Jacobian。 */
  Matrix<Scalar, FDim, XDim> dfdx;
  /** @brief 对输入的 Jacobian。 */
  Matrix<Scalar, FDim, UDim> dfdu;
  /** @brief 常数项向量。 */
  Vector<Scalar, FDim> f;

  /** @brief 默认构造。 */
  VectorFunctionQuadraticApproximation() = default;

  /** @brief 将各系数置零。 */
  VectorFunctionQuadraticApproximation& setZero() {
    for (int i = 0; i < FDim; ++i) {
      dfdxx[i].setZero();
      dfdux[i].setZero();
      dfduu[i].setZero();
    }
    dfdx.setZero();
    dfdu.setZero();
    f.setZero();
    return *this;
  }

  /** @brief 返回零初始化的近似对象。 */
  static VectorFunctionQuadraticApproximation Zero() {
    VectorFunctionQuadraticApproximation approximation;
    approximation.setZero();
    return approximation;
  }
};

template <typename Scalar, int FDim, int XDim>
struct VectorFunctionQuadraticApproximation<Scalar, FDim, XDim, 0> {
  /** @brief 每个输出分量对状态的 Hessian。 */
  std::array<Matrix<Scalar, XDim, XDim>, FDim> dfdxx;
  /** @brief 对状态的 Jacobian。 */
  Matrix<Scalar, FDim, XDim> dfdx;
  /** @brief 常数项向量。 */
  Vector<Scalar, FDim> f;

  /** @brief 默认构造。 */
  VectorFunctionQuadraticApproximation() = default;

  /** @brief 将各系数置零。 */
  VectorFunctionQuadraticApproximation& setZero() {
    for (int i = 0; i < FDim; ++i) {
      dfdxx[i].setZero();
    }
    dfdx.setZero();
    f.setZero();
    return *this;
  }

  /** @brief 返回零初始化的近似对象。 */
  static VectorFunctionQuadraticApproximation Zero() {
    VectorFunctionQuadraticApproximation approximation;
    approximation.setZero();
    return approximation;
  }
};