/**
 * @file Initializer.hpp
 * @brief 轨迹初始化器接口：在无控制器的时间段内，由 (time, state) 与 nextTime
 * 计算 input 与 nextState。
 */
#pragma once

#include "Matrix/Types.hpp"

/**
 * @brief 求解器在无控制器可用时使用的初始化器接口；简单实现见
 * DefaultInitializer。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class Initializer {
 public:
  Initializer() = default;
  virtual ~Initializer() = default;

  /**
   * @brief 由当前时间与状态及下一时刻，计算当前输入与下一状态。
   * @param [in] time 当前时间。
   * @param [in] state 当前状态。
   * @param [in] nextTime 下一时刻（通常为 time + timeStep）。
   * @param [out] input 当前段的输入。
   * @param [out] nextState 下一时刻的状态。
   */
  virtual void compute(const Scalar time, const Vector<Scalar, XDim>& state,
                       const Scalar nextTime, Vector<Scalar, UDim>& input,
                       Vector<Scalar, XDim>& nextState) = 0;

 protected:
  /** @brief 拷贝构造（保护）。 */
  Initializer(const Initializer& rhs) = default;
};
