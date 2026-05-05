/**
 * @file PrimalSolution.hpp
 * @brief 原始问题解：时间/状态/输入轨迹及线性控制器。
 */
#pragma once
#include <array>
#include "LinearController.hpp"

/**
 * @brief 原始问题解：一条 rollout 的时间、状态、输入轨迹及对应的线性控制器。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 * @tparam PredictLength 预测步数（轨迹点数为 PredictLength+1）。
 */
template <typename Scalar, int XDim, int UDim, size_t PredictLength>
struct PrimalSolution
{
  /** @brief 默认构造。 */
  PrimalSolution() = default;

  /** @brief 析构函数。 */
  ~PrimalSolution() = default;

  /** @brief 拷贝构造。 */
  PrimalSolution(const PrimalSolution &other)
      : timeTrajectory_(other.timeTrajectory_),
        stateTrajectory_(other.stateTrajectory_),
        inputTrajectory_(other.inputTrajectory_),
        controller_(other.controller_)
  {
  }

  /** @brief 拷贝赋值。 */
  PrimalSolution &operator=(const PrimalSolution &other)
  {
    timeTrajectory_ = other.timeTrajectory_;
    stateTrajectory_ = other.stateTrajectory_;
    inputTrajectory_ = other.inputTrajectory_;
    controller_ = other.controller_;
    return *this;
  }

  /** @brief 移动构造。 */
  PrimalSolution(PrimalSolution &&other) noexcept = default;

  /** @brief 移动赋值。 */
  PrimalSolution &operator=(PrimalSolution &&other) noexcept = default;

  /** @brief 与另一 PrimalSolution 交换时间/状态/输入轨迹及控制器。 */
  void swap(PrimalSolution &other)
  {
    timeTrajectory_.swap(other.timeTrajectory_);
    stateTrajectory_.swap(other.stateTrajectory_);
    inputTrajectory_.swap(other.inputTrajectory_);
    controller_.swap(other.controller_);
  }

  /** @brief 清空：控制器清空，时间/状态/输入轨迹置零。 */
  void clear()
  {
    controller_.clear();
    for (size_t i = 0; i < PredictLength + 1; ++i)
    {
      timeTrajectory_[i] = 0;
      stateTrajectory_[i].setZero();
      inputTrajectory_[i].setZero();
    }
  }

  /** @brief 时间序列，长度 PredictLength+1。 */
  std::array<Scalar, PredictLength + 1> timeTrajectory_;
  /** @brief 状态轨迹。 */
  std::array<Vector<Scalar, XDim>, PredictLength + 1> stateTrajectory_;
  /** @brief 输入轨迹。 */
  std::array<Vector<Scalar, UDim>, PredictLength + 1> inputTrajectory_;
  /** @brief 线性控制器（时间戳与增益/偏置数组）。 */
  LinearController<Scalar, XDim, UDim, PredictLength + 1> controller_;
};
