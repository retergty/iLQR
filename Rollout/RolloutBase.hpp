/**
 * @file RolloutBase.hpp
 * @brief 前向 rollout 基类：在 [initTime, finalTime]
 * 上按给定控制器展开状态轨迹并写入缓冲区。
 */
#pragma once

#include "Controller/Controller.hpp"
#include "Types.hpp"

/** @brief 求根算法类型枚举（当前未使用）。 */
enum class RootFinderType { ANDERSON_BJORCK, PEGASUS, ILLINOIS, REGULA_FALSI };

/**
 * @brief 前向 rollout 的配置：时间步长、是否重建输入轨迹。
 */
template <typename Scalar>
struct RolloutSettings {
  /** @brief 固定步长 rollout 使用的时间步长。 */
  Scalar timeStep = 1e-2;

  /** @brief rollout 时是否写入或重建输入轨迹。 */
  bool reconstructInputTrajectory = true;
};

/**
 * @brief 指向轨迹缓冲区的轻量句柄：时间/状态/输入指针及最大长度。
 */
template <typename Scalar, int XDim, int UDim>
struct RolloutTrajectoryPointer {
  /** @brief 构造：绑定时间、状态、输入数组指针及最大写入长度。 */
  RolloutTrajectoryPointer(Scalar* time_trajectory,
                           Vector<Scalar, XDim>* state_trajectory,
                           Vector<Scalar, UDim>* input_trajectory,
                           int max_length)
      : timeTrajectory(time_trajectory),
        stateTrajectory(state_trajectory),
        inputTrajectory(input_trajectory),
        maxLength(max_length) {};
  Scalar* timeTrajectory;
  Vector<Scalar, XDim>* stateTrajectory;
  Vector<Scalar, UDim>* inputTrajectory;
  size_t maxLength;
};

/**
 * @brief 前向 rollout 抽象基类：用给定控制器与初态在 [initTime, finalTime]
 * 上生成状态、输入和时间轨迹。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class RolloutBase {
 public:
  using RolloutTrajectoryPointer_t =
      RolloutTrajectoryPointer<Scalar, XDim, UDim>;
  /** @brief 默认构造。 */
  explicit RolloutBase() {}

  /** @brief 虚析构。 */
  virtual ~RolloutBase() = default;

  /** @brief 返回 rollout 配置（步长、是否重建输入轨迹等）。 */
  RolloutSettings<Scalar>& settings() { return rolloutSettings_; }

  /**
   * @brief 从 initTime 到 finalTime 用给定控制器前向展开轨迹。
   * @param [in] initTime 初始时间。
   * @param [in] initState 初始状态。
   * @param [in] finalTime 终止时间。
   * @param [in] controller 控制策略（可为空）。
   * @param [in,out] trajectory 轨迹输出（时间/状态/输入指针及最大长度）。
   * @return 写入的轨迹点数。
   */
  virtual int run(const Scalar initTime, const Vector<Scalar, XDim>& initState,
                  const Scalar finalTime,
                  ControllerBase<Scalar, XDim, UDim>* controller,
                  RolloutTrajectoryPointer_t& trajectory) = 0;

 protected:
  RolloutSettings<Scalar> rolloutSettings_{};
};