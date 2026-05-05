/**
 * @file Reference.hpp
 * @brief 参考轨迹接口：用户定义的目标时间/状态/输入轨迹及插值查询。
 */
#pragma once
#include "Types.hpp"
#include "Integration.hpp"
#include "LinearInterpolation.hpp"
#include <array>
#include <algorithm>

/**
 * @brief 目标轨迹容器：时间、期望状态与期望输入的离散序列，支持按索引或时间插值查询。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDim, int UDim, int ArrayLength>
struct TargetTrajectories
{
  /** @brief 默认构造。 */
  TargetTrajectories() = default;

  /** @brief 用给定时间、状态、输入轨迹构造。 */
  TargetTrajectories(
    const std::array<Scalar, ArrayLength>& desiredTimeTrajectory,
    const std::array<Vector<Scalar, XDim>, ArrayLength>& desiredStateTrajectory,
    const std::array<Vector<Scalar, UDim>, ArrayLength>& desiredInputTrajectory)
    : timeTrajectory(desiredTimeTrajectory), stateTrajectory(desiredStateTrajectory), inputTrajectory(desiredInputTrajectory) {
  }

  /** @brief 用给定时间、状态轨迹构造，输入轨迹置零。 */
  TargetTrajectories(
    const std::array<Scalar, ArrayLength>& desiredTimeTrajectory,
    const std::array<Vector<Scalar, XDim>, ArrayLength>& desiredStateTrajectory)
    : timeTrajectory(desiredTimeTrajectory), stateTrajectory(desiredStateTrajectory)
  {
    for (auto& vec : inputTrajectory)
    {
      vec.setZero();
    }
  }

  /** @brief 设置完整轨迹（时间、状态、输入）。 */
  void setTrajectory(const std::array<Scalar, ArrayLength>& desiredTimeTrajectory,
    const std::array<Vector<Scalar, XDim>, ArrayLength>& desiredStateTrajectory,
    const std::array<Vector<Scalar, UDim>, ArrayLength>& desiredInputTrajectory)
  {
    timeTrajectory = desiredTimeTrajectory;
    stateTrajectory = desiredStateTrajectory;
    inputTrajectory = desiredInputTrajectory;
  }

  /** @brief 设置时间与状态轨迹，输入轨迹置零。 */
  void setTrajectory(const std::array<Scalar, ArrayLength>& desiredTimeTrajectory,
    const std::array<Vector<Scalar, XDim>, ArrayLength>& desiredStateTrajectory)
  {
    timeTrajectory = desiredTimeTrajectory;
    stateTrajectory = desiredStateTrajectory;
    for (auto& vec : inputTrajectory)
    {
      vec.setZero();
    }
  }

  /** @brief 按索引获取期望状态。 */
  Vector<Scalar, XDim> getDesiredState(const int index) const
  {
    assert(index >= 0 && index < ArrayLength);
    return stateTrajectory[index];
  }
  /** @brief 按编译期索引获取期望状态。 */
  template <int Index>
  Vector<Scalar, XDim> getDesiredState() const
  {
    assert(Index >= 0 && Index < ArrayLength);
    return stateTrajectory[Index];
  }
  /** @brief 按时间插值获取期望状态。 */
  Vector<Scalar, XDim> getDesiredState(const Scalar time) const
  {
    return LinearInterpolation::interpolate(time, timeTrajectory, stateTrajectory);
  }

  /** @brief 按索引获取期望输入。 */
  Vector<Scalar, UDim> getDesiredInput(const int index) const
  {
    assert(index >= 0 && index < ArrayLength);
    return inputTrajectory[index];
  }

  /** @brief 按编译期索引获取期望输入。 */
  template <int Index>
  Vector<Scalar, UDim> getDesiredInput() const
  {
    assert(Index >= 0 && Index < ArrayLength);
    return inputTrajectory[Index];
  }
  /** @brief 按时间插值获取期望输入。 */
  Vector<Scalar, UDim> getDesiredInput(const Scalar time) const
  {
    return LinearInterpolation::interpolate(time, timeTrajectory, inputTrajectory);
  }

  /** @brief 时间序列。 */
  std::array<Scalar, ArrayLength> timeTrajectory;
  /** @brief 期望状态轨迹。 */
  std::array<Vector<Scalar, XDim>, ArrayLength> stateTrajectory;
  /** @brief 期望输入轨迹。 */
  std::array<Vector<Scalar, UDim>, ArrayLength> inputTrajectory;
};