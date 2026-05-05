/**
 * @file SystemDynamicsBase.hpp
 * @brief 系统动力学与线性化基类：定义线性化流映射与跳变映射。
 *
 * 线性化流映射：\f$ dx/dt = A(t) \delta x + B(t) \delta u \f$
 * 线性化跳变映射：\f$ x^+ = G \delta x + H \delta u \f$
 */
#pragma once
#include "ControlledSystemBase.hpp"

/**
 * @brief 系统动力学与线性化基类：子类需实现 linearApproximation。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam UDimisions 控制维度。
 */
template <typename Scalar, int XDimisions, int UDimisions>
class SystemDynamicsBase : public ControlledSystemBase<Scalar, XDimisions, UDimisions>
{
public:
  /** @brief 默认构造。 */
  SystemDynamicsBase() = default;

  /** @brief 虚析构函数。 */
  ~SystemDynamicsBase() override = default;

  /**
   * @brief 计算动力学的线性近似（雅可比 dfdx、dfdu 与 f）。
   * @param [in] t 当前时间。
   * @param [in] x 当前状态。
   * @param [in] u 当前输入。
   * @return 状态对时间导数的线性近似（f, dfdx, dfdu）。
   */
  virtual VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>
  linearApproximation(Scalar t, const Vector<Scalar, XDimisions> &x, const Vector<Scalar, UDimisions> &u) = 0;
};
