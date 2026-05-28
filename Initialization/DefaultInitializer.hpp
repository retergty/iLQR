/**
 * @file DefaultInitializer.hpp
 * @brief 默认初始化器：输入置零，下一状态等于当前状态（保持不动）。
 */
#pragma once

#include "Initialization/Initializer.hpp"

/**
 * @brief 默认初始化器实现：input 置零，nextState = state。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class DefaultInitializer final : public Initializer<Scalar, XDim, UDim> {
 public:
  /** @brief 默认构造。 */
  explicit DefaultInitializer() = default;

  ~DefaultInitializer() override = default;

  /** @brief 将 input 置零，nextState 设为当前 state。 */
  void compute(const Scalar time, const Vector<Scalar, XDim>& state,
               const Scalar nextTime, Vector<Scalar, UDim>& input,
               Vector<Scalar, XDim>& nextState) override {
    (void)time;
    (void)nextTime;
    (void)state;

    input.setZero();
    nextState = state;
  }

 protected:
  DefaultInitializer(const DefaultInitializer& rhs) = default;
};
