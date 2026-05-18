/**
 * @file TimeTriggeredRollout.hpp
 * @brief 按固定时间步长做前向 rollout：用 RK45
 * 积分受控动力学并可选地重建输入轨迹。
 */
#pragma once
#include "ControlledSystemBase.hpp"
#include "RolloutBase.hpp"
#include "RungeKuttaDormandPrince5.hpp"

/**
 * @brief 时间触发前向 rollout：绑定受控系统动力学，用 RK45
 * 固定步长积分，结果写入 trajectory。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class TimeTriggeredRollout : public RolloutBase<Scalar, XDim, UDim> {
 public:
  using RolloutTrajectoryPointer_t =
      typename RolloutBase<Scalar, XDim, UDim>::RolloutTrajectoryPointer_t;
  /**
   * @brief 构造：绑定系统动力学与积分步长。
   * @param [in] systemDynamics 用于前向积分的受控系统动力学指针。
   * @param [in] timeStep 积分步长。
   */
  explicit TimeTriggeredRollout(
      ControlledSystemBase<Scalar, XDim, UDim>* systemDynamics,
      const Scalar timeStep)
      : systemDynamicsPtr_(systemDynamics) {
    this->rolloutSettings_.timeStep = timeStep;
  };

  ~TimeTriggeredRollout() override = default;

  /** @brief 返回底层动力学指针。 */
  ControlledSystemBase<Scalar, XDim, UDim>* systemDynamicsPtr() {
    return systemDynamicsPtr_;
  }

  /**
   * @brief 用当前控制器从 initTime 积分到 finalTime，状态/时间写入
   * trajectory，可选再算输入轨迹。
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

    // set controller
    systemDynamicsPtr_->setController(controller);

    Observer<Scalar, XDim> observer(
        trajectory.maxLength, trajectory.stateTrajectory,
        trajectory.timeTrajectory);  // concatenate trajectory
    // integrate controlled system
    RK45Intergraor_.integrateConst(*systemDynamicsPtr_, observer, initState,
                                   initTime, finalTime,
                                   this->settings().timeStep);

    int RolloutIntegrateCount = observer.getCount();

    // compute control input trajectory and concatenate to inputTrajectory
    if (this->settings().reconstructInputTrajectory) {
      for (int i = 0; i < RolloutIntegrateCount; i++) {
        trajectory.inputTrajectory[i] =
            systemDynamicsPtr_->controllerPtr()->computeInput(
                trajectory.timeTrajectory[i], trajectory.stateTrajectory[i]);
      }  // end of k_u loop
    }

    return RolloutIntegrateCount;
  }

  // Scalar run(const Scalar initTime, const Vector<Scalar, XDim>& initState,
  // const int steps, ControllerBase<Scalar, XDim, UDim>& controller,
  //     std::array<Scalar, ArrayLen>& timeTrajectory, std::array<Vector<Scalar,
  //     XDim>, ArrayLen>& stateTrajectory, std::array<Vector<Scalar, UDim>,
  //     ArrayLen>& inputTrajectory) override
  // {
  //     Scalar finalTime = initTime + steps * this->settings().timeStep;
  //     // set controller
  //     systemDynamicsPtr_->setController(&controller);

  //     Observer<Scalar, XDim> observer(ArrayLen, stateTrajectory.data(),
  //     timeTrajectory.data()); // concatenate trajectory
  //     // integrate controlled system
  //     RK45Intergraor_.integrateConst(*systemDynamicsPtr_, observer,
  //     initState, initTime, finalTime, this->settings().timeStep);
  // }

 private:
  ControlledSystemBase<Scalar, XDim, UDim>* systemDynamicsPtr_{nullptr};

  RungeKuttaDormandPrince5<Scalar, XDim> RK45Intergraor_;
};
