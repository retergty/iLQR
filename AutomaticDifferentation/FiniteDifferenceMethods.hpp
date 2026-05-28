/**
 * @file FiniteDifferenceMethods.hpp
 * @brief 有限差分数值求导工具：用于计算向量函数对状态或输入的 Jacobian。
 */
#pragma once
#include <math.h>

#include <algorithm>

#include "Dynamics/ControlledSystemBase.hpp"
#include "Types.hpp"

/**
 * @brief 用有限差分计算向量函数对变量的 Jacobian。
 *
 * 步长会按变量当前量级缩放，以降低固定步长在大/小数值下的舍入误差。
 * 当 `doubleSidedDerivative=true` 时使用中心差分，否则使用前向差分。
 *
 * @tparam Scalar 标量类型。
 * @tparam StateDimisions 函数输出维度。
 * @tparam VarDimisions 自变量维度。
 * @tparam Function 可调用对象类型，签名应兼容 `Vector -> Vector`。
 * @param [in] f 待求导函数。
 * @param [in] x0 线性化点。
 * @param [in] eps 基础扰动尺度。
 * @param [in] doubleSidedDerivative 是否使用中心差分。
 * @return 在 `x0` 处计算得到的 Jacobian。
 */
template <typename Scalar, int StateDimisions, int VarDimisions,
          typename Function>
Matrix<Scalar, StateDimisions, VarDimisions> finiteDifferenceDerivative(
    const Function& f, const Vector<Scalar, VarDimisions>& x0, Scalar eps,
    bool doubleSidedDerivative) {
  const Vector<Scalar, StateDimisions> f0 = f(x0);
  Matrix<Scalar, StateDimisions, VarDimisions> jacobian;

  for (size_t i = 0; i < VarDimisions; i++) {
    // 灵感来源：
    // http://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating_point_arithmetic
    Scalar h = eps * std::max(std::fabs(x0(i)), 1.0);

    Vector<Scalar, VarDimisions> xPlus = x0;
    xPlus(i) += h;

    if (doubleSidedDerivative) {
      Vector<Scalar, VarDimisions> xMinus = x0;
      xMinus(i) -= h;
      jacobian.col(i) = (f(xPlus) - f(xMinus)) / (2.0 * h);
    } else {
      jacobian.col(i) = (f(xPlus) - f0) / h;
    }
  }

  return jacobian;
}

/**
 * @brief 用有限差分计算动力学对状态的 Jacobian。
 *
 * 当 `isSecondOrderSystem=true` 时，会将上半部分修正为二阶系统常见的
 * `d/dx [q_dot] = 0` 与 `d/dv [q_dot] = I` 结构，假设状态排列为 `[q, q_dot]`。
 *
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @param [in] system 连续时间受控系统。
 * @param [in] t 当前时间。
 * @param [in] x 线性化状态。
 * @param [in] u 线性化输入。
 * @param [in] eps 基础扰动尺度。
 * @param [in] doubleSidedDerivative 是否使用中心差分。
 * @param [in] isSecondOrderSystem 是否按二阶系统结构修正 Jacobian。
 * @return 状态 Jacobian `A = d f / d x`。
 */
template <typename Scalar, int XDim, int UDim>
Matrix<Scalar, XDim, XDim> finiteDifferenceDerivativeState(
    ControlledSystemBase<Scalar, XDim, UDim>& system, Scalar t,
    const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u, Scalar eps,
    bool doubleSidedDerivative, bool isSecondOrderSystem) {
  auto f = [&](const Vector<Scalar, XDim>& var) -> Vector<Scalar, XDim> {
    return system.computeFlowMap(t, var, u);
  };

  Matrix<Scalar, XDim, XDim> A = finiteDifferenceDerivative<Scalar, XDim, XDim>(
      f, x, eps, doubleSidedDerivative);

  if (isSecondOrderSystem) {
    // 假定状态向量 = [x, x_dot]。
    constexpr int HalfXDim = XDim / 2;
    A.template topLeftCorner<HalfXDim, HalfXDim>().setZero();
    A.template topRightCorner<HalfXDim, HalfXDim>().setIdentity();
  }
  return A;
}

/**
 * @brief 用有限差分计算动力学对输入的 Jacobian。
 *
 * 当 `isSecondOrderSystem=true` 时，会将上半部分置零，反映
 * `[q, q_dot]` 形式下位置导数通常不直接依赖控制输入的结构假设。
 *
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @param [in] system 连续时间受控系统。
 * @param [in] t 当前时间。
 * @param [in] x 线性化状态。
 * @param [in] u 线性化输入。
 * @param [in] eps 基础扰动尺度。
 * @param [in] doubleSidedDerivative 是否使用中心差分。
 * @param [in] isSecondOrderSystem 是否按二阶系统结构修正 Jacobian。
 * @return 输入 Jacobian `B = d f / d u`。
 */
template <typename Scalar, int XDim, int UDim>
Matrix<Scalar, XDim, UDim> finiteDifferenceDerivativeInput(
    ControlledSystemBase<Scalar, XDim, UDim>& system, Scalar t,
    const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u, Scalar eps,
    bool doubleSidedDerivative, bool isSecondOrderSystem) {
  auto f = [&](const Vector<Scalar, UDim>& var) -> Vector<Scalar, XDim> {
    return system.computeFlowMap(t, x, var);
  };

  Matrix<Scalar, XDim, UDim> B = finiteDifferenceDerivative<Scalar, XDim, UDim>(
      f, u, eps, doubleSidedDerivative);

  if (isSecondOrderSystem) {
    // 假定状态向量 = [x, x_dot]。
    constexpr int HalfXDim = XDim / 2;
    B.template topRows<HalfXDim>().setZero();
  }
  return B;
}