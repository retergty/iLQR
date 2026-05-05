/**
 * @file ControlledSystemBase.hpp
 * @brief 受控系统动力学基类：支持按 (t,x) 或 (t,x,u) 计算状态导数，并绑定控制器。
 */
#pragma once
#include "Types.hpp"
#include "Integration.hpp"
#include "Controller.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"

/**
 * @brief 受控系统动力学基类：可设置控制器，computeFlowMap(t,x) 内部用控制器算 u 再调用 computeFlowMap(t,x,u)。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam UDimisions 控制维度。
 */
template <typename Scalar, int XDimisions, int UDimisions>
class ControlledSystemBase : public OdeBase<Scalar, XDimisions>
{
public:
  /** @brief 默认构造。 */
  ControlledSystemBase() = default;

  /** @brief 虚析构。 */
  ~ControlledSystemBase() override = default;

  /**
   * @brief 用当前控制器根据 (t,x) 计算 u，再计算状态导数 dx/dt。
   * @param [in] t 当前时间。
   * @param [in] x 当前状态。
   * @return 状态对时间的导数。
   */
  Vector<Scalar, XDimisions> computeFlowMap(Scalar t, const Vector<Scalar, XDimisions> &x) const override final
  {
    assert(controllerPtr_ != nullptr);
    const Vector<Scalar, UDimisions> u = controllerPtr_->computeInput(t, x);
    return computeFlowMap(t, x, u);
  }

  /**
   * @brief 给定 (t,x,u) 计算状态导数 dx/dt（子类实现）。
   * @param [in] t 当前时间。
   * @param [in] x 当前状态。
   * @param [in] u 当前输入。
   * @return 状态对时间的导数。
   */
  virtual Vector<Scalar, XDimisions> computeFlowMap(Scalar t, const Vector<Scalar, XDimisions> &x, const Vector<Scalar, UDimisions> &u) const = 0;

  /** @brief 设置控制器指针，用于 rollout 时计算 u。 */
  void setController(ControllerBase<Scalar, XDimisions, UDimisions> *controllerPtr)
  {
    controllerPtr_ = controllerPtr;
  }
  /** @brief 返回当前绑定的控制器指针。 */
  const ControllerBase<Scalar, XDimisions, UDimisions> *controllerPtr() const
  {
    return controllerPtr_;
  }

private:
  ControllerBase<Scalar, XDimisions, UDimisions> *controllerPtr_ = nullptr;
};
