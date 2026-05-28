/**
 * @file InitializerRollout.hpp
 * @brief 用 Initializer 生成轨迹的 rollout：不积分动力学，按步长调用
 * initializer.compute 得到状态与输入。
 */
#pragma once

#include "Initialization/Initializer.hpp"
#include "Rollout/RolloutBase.hpp"

/**
 * @brief 基于初始化器的 rollout：在 [initTime, finalTime] 上按固定步长调用
 * initializer，填充时间/状态/输入轨迹。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class InitializerRollout : RolloutBase<Scalar, XDim, UDim> {
 public:
  using Initializer_t = Initializer<Scalar, XDim, UDim>;
  using RolloutTrajectoryPointer_t =
      typename RolloutBase<Scalar, XDim, UDim>::RolloutTrajectoryPointer_t;
  /**
   * @brief 构造：绑定初始化器与步长。
   * @param [in] initializer 用于生成状态与输入的初始化器。
   * @param [in] timeStep 时间步长。
   */
  explicit InitializerRollout(Initializer_t& initializer, const Scalar timeStep)
      : initializer_(initializer) {
    this->rolloutSettings_.timeStep = timeStep;
  }

  ~InitializerRollout() override = default;

  /**
   * @brief 从 initTime 到 finalTime 按步长调用
   * initializer.compute，将时间/状态/输入写入 trajectory。
   * @param [in] initTime 初始时间。
   * @param [in] initState 初始状态。
   * @param [in] finalTime 终止时间。
   * @param [in] controller 未使用。
   * @param [in,out] trajectory 输出轨迹。
   * @return 步数（写入点数减 1）。
   */
  int run(const Scalar initTime, const Vector<Scalar, XDim>& initState,
          const Scalar finalTime,
          ControllerBase<Scalar, XDim, UDim>* controller,
          RolloutTrajectoryPointer_t& trajectory) override {
    assert(finalTime > initTime);
    (void)controller;

    // 通过加入一个 dt 的小分数确保包含 finalTime，使得：N
    // * dt <= finalTime < (N + 1) * dt。
    Scalar finalTimeLocal = finalTime + 0.1 * this->settings().timeStep;
    Scalar t = initTime;
    const Scalar timeStep = this->settings().timeStep;

    Vector<Scalar, XDim> state = initState;
    Vector<Scalar, XDim> nextState;

    const size_t numSteps = (finalTimeLocal - initTime) / timeStep;

    for (size_t i = 0; i < numSteps + 1; i++) {
      assert(i < trajectory.maxLength);
      initializer_.compute(t, state, t + timeStep,
                           trajectory.inputTrajectory[i], nextState);
      trajectory.timeTrajectory[i] = t;
      trajectory.stateTrajectory[i] = state;
      state = nextState;
      t += timeStep;
    }  // i 循环结束。
    return numSteps;
  }

 private:
  Initializer_t& initializer_;
};