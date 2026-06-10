/**
 * @file DiscreteTimeRollout.hpp
 * @brief 按固定时间步长做离散系统前向 rollout：逐步调用离散状态转移映射。
 */
#pragma once

#include <cassert>

#include "Controller/Controller.hpp"
#include "Dynamics/DiscreteSystemBase.hpp"
#include "Rollout/RolloutBase.hpp"

/**
 * @brief 离散时间前向 rollout：绑定离散系统动力学，用 computeMap()
 * 逐节点推进状态，结果写入 trajectory。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class DiscreteTimeRollout : public RolloutBase<Scalar, XDim, UDim> {
 public:
  using RolloutTrajectoryPointer_t =
      typename RolloutBase<Scalar, XDim, UDim>::RolloutTrajectoryPointer_t;

  /**
   * @brief 构造：绑定离散系统与采样步长。
   * @param [in] systemDynamics 用于前向 rollout 的离散系统指针。
   * @param [in] timeStep 采样步长。
   */
  explicit DiscreteTimeRollout(
      DiscreteSystemBase<Scalar, XDim, UDim>* systemDynamics,
      const Scalar timeStep)
      : systemDynamicsPtr_(systemDynamics) {
    this->rolloutSettings_.timeStep = timeStep;
  }

  ~DiscreteTimeRollout() override = default;

  /** @brief 返回底层离散系统指针。 */
  DiscreteSystemBase<Scalar, XDim, UDim>* systemDynamicsPtr() {
    return systemDynamicsPtr_;
  }

  /**
   * @brief 用当前控制器从 initTime 到 finalTime 逐步推进离散系统。
   * @param [in] initTime 初始时间。
   * @param [in] initState 初始状态。
   * @param [in] finalTime 终止时间。
   * @param [in] controller 控制策略。
   * @param [in,out] trajectory 输出轨迹缓冲区。
   * @return 写入的轨迹点数。
   */
  int run(const Scalar initTime, const Vector<Scalar, XDim>& initState,
          const Scalar finalTime,
          ControllerBase<Scalar, XDim, UDim>* controller,
          RolloutTrajectoryPointer_t& trajectory) override {
    assert(finalTime > initTime);
    assert(systemDynamicsPtr_ != nullptr);
    assert(controller != nullptr);

    const auto& settings = this->settings();
    const Scalar timeStep = settings.timeStep;
    const bool reconstructInputTrajectory = settings.reconstructInputTrajectory;
    const Scalar finalTimeLocal = finalTime + Scalar(0.1) * timeStep;
    const size_t numSteps =
        static_cast<size_t>((finalTimeLocal - initTime) / timeStep);
    assert(numSteps < trajectory.maxLength);

    Scalar t = initTime;
    trajectory.timeTrajectory[0] = t;
    trajectory.stateTrajectory[0] = initState;

    for (size_t i = 0; i < numSteps; ++i) {
      const Vector<Scalar, UDim> input =
          controller->computeInput(t, trajectory.stateTrajectory[i]);

      if (reconstructInputTrajectory) {
        trajectory.inputTrajectory[i] = input;
      }

      trajectory.stateTrajectory[i + 1] = systemDynamicsPtr_->computeMap(
          t, trajectory.stateTrajectory[i], input, timeStep);
      t += timeStep;
      trajectory.timeTrajectory[i + 1] = t;
    }

    if (reconstructInputTrajectory) {
      trajectory.inputTrajectory[numSteps] =
          controller->computeInput(t, trajectory.stateTrajectory[numSteps]);
    }

    return static_cast<int>(numSteps + 1);
  }

 private:
  DiscreteSystemBase<Scalar, XDim, UDim>* systemDynamicsPtr_{nullptr};
};
