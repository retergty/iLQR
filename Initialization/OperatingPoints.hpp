/**
 * @file OperatingPoints.hpp
 * @brief 工作点初始化器：用固定状态/输入工作点生成轨迹，不依赖动力学积分。
 */
#pragma once

#include "Initializer.hpp"

/**
 * @brief 基于工作点的初始化器：输出恒为给定状态工作点与输入工作点，nextState
 * 为状态工作点。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class OperatingPoints final : public Initializer<Scalar, XDim, UDim> {
 public:
  /**
   * @brief 用状态工作点与输入工作点构造。
   * @param [in] stateOperatingPoint 状态工作点。
   * @param [in] inputOperatingPoint 输入工作点。
   */
  OperatingPoints(const Vector<Scalar, XDim>& stateOperatingPoint,
                  const Vector<Scalar, UDim>& inputOperatingPoint)
      : stateTrajectory_(stateOperatingPoint),
        inputTrajectory_(inputOperatingPoint) {}

  /** @brief 析构函数。 */
  ~OperatingPoints() override = default;

  /** @brief 将 input 设为输入工作点，nextState 设为状态工作点。 */
  void compute(const Scalar time, const Vector<Scalar, XDim>& state,
               const Scalar nextTime, Vector<Scalar, UDim>& input,
               Vector<Scalar, XDim>& nextState) override {
    (void)time;
    (void)state;
    (void)nextTime;
    input = inputTrajectory_;
    nextState = stateTrajectory_;
  }

 private:
  /** @brief 拷贝构造（保护）。 */
  OperatingPoints(const OperatingPoints& other) = default;

  /** @brief 状态工作点。 */
  const Vector<Scalar, XDim> stateTrajectory_;
  /** @brief 输入工作点。 */
  const Vector<Scalar, UDim> inputTrajectory_;
};
