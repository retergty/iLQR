/**
 * @file Integration.hpp
 * @brief 常微分方程积分接口：自治系统基类 OdeBase、积分器类型与
 * IntegratorBase。
 */
#pragma once
#include "Observer.hpp"
#include "Types.hpp"

/**
 * @brief 自治系统动力学基类：根据时间与状态计算状态导数 dx/dt。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 */
template <typename Scalar, int XDim>
class OdeBase {
 public:
  /** @brief 默认构造。 */
  OdeBase() = default;

  /** @brief 虚析构。 */
  virtual ~OdeBase() = default;

  /**
   * @brief 计算当前时刻的状态导数。
   * @param [in] t 当前时间。
   * @param [in] x 当前状态。
   * @return 状态对时间的导数。
   */
  virtual Vector<Scalar, XDim> computeFlowMap(
      Scalar t, const Vector<Scalar, XDim>& x) = 0;
};

/**
 * @brief 积分器类型枚举：欧拉、ODE45、RK4 等。
 */
enum class IntegratorType {
  EULER,
  ODE45,
  ODE45_OCS2,
  ADAMS_BASHFORTH,
  BULIRSCH_STOER,
  MODIFIED_MIDPOINT,
  RK4,
  RK5_VARIABLE,
  ADAMS_BASHFORTH_MOULTON
};

/**
 * @brief 自治系统积分器基类：固定步长积分，将轨迹写入 Observer。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 */
template <typename Scalar, int XDim>
class IntegratorBase {
 public:
  /** @brief 默认构造。 */
  IntegratorBase() = default;

  /** @brief 虚析构。 */
  virtual ~IntegratorBase() = default;

  /**
   * @brief 从 startTime 到 finalTime 以固定步长 dt 积分，结果通过 observer
   * 输出。
   * @param [in] system 系统动力学。
   * @param [in,out] observer 观测器（接收时间与状态）。
   * @param [in] initialState 初始状态。
   * @param [in] startTime 初始时间。
   * @param [in] finalTime 终止时间。
   * @param [in] dt 时间步长。
   */
  virtual void integrateConst(OdeBase<Scalar, XDim>& system,
                              Observer<Scalar, XDim>& observer,
                              const Vector<Scalar, XDim>& initialState,
                              const Scalar startTime, const Scalar finalTime,
                              const Scalar dt) = 0;
};