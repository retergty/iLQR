/**
 * @file SystemDynamicsLinearizer.hpp
 * @brief 系统动力学数值线性化器：通过有限差分构造 `f`, `dfdx`, `dfdu`。
 */
#pragma once

#include <memory>

#include "ControlledSystemBase.hpp"
#include "SystemDynamicsBase.hpp"
#include "FiniteDifferenceMethods.hpp"

/**
 * @brief 系统动力学数值线性化器。
 *
 * 给定一个非线性受控系统，在指定 `(t, x, u)` 处通过有限差分构造
 * 连续时间线性近似：
 *
 * \f[
 * \dot{x} \approx f(t, x, u) + A(t)\delta x + B(t)\delta u
 * \f]
 *
 * 其中 `A = d f / d x`，`B = d f / d u`。
 *
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam UDimisions 输入维度。
 */
template <typename Scalar, int XDimisions, int UDimisions>
class SystemDynamicsLinearizer final : public SystemDynamicsBase<Scalar, XDimisions, UDimisions> {
public:
  /** @brief 状态向量类型。 */
  using StateVector_t = Vector<Scalar, XDimisions>;
  /** @brief 输入向量类型。 */
  using InputVector_t = Vector<Scalar, UDimisions>;
  /** @brief 动力学线性近似类型。 */
  using VectorFunctionLinearApproximation_t = VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>;

  /**
   * @brief 构造数值线性化器。
   * @param [in] nonlinearSystemPtr 被线性化的系统指针，生命周期由外部管理。
   * @param [in] doubleSidedDerivative 是否使用中心差分。
   * @param [in] isSecondOrderSystem 是否按二阶系统结构修正 Jacobian。
   * @param [in] eps 基础有限差分步长。
   */
  explicit SystemDynamicsLinearizer(ControlledSystemBase<Scalar, XDimisions, UDimisions>* nonlinearSystemPtr, bool doubleSidedDerivative = true,
    bool isSecondOrderSystem = false, Scalar eps = Eigen::NumTraits<Scalar>::epsilon()) : SystemDynamicsBase<Scalar, XDimisions, UDimisions>(),
    controlledSystemPtr_(nonlinearSystemPtr),
    doubleSidedDerivative_(doubleSidedDerivative),
    isSecondOrderSystem_(isSecondOrderSystem),
    eps_(eps)
  {
  }

  ~SystemDynamicsLinearizer() override = default;

  /** @brief 直接转发到底层非线性系统的流映射。 */
  StateVector_t computeFlowMap(Scalar time, const StateVector_t& state, const InputVector_t& input)const override
  {
    return controlledSystemPtr_->computeFlowMap(time, state, input);
  }

  /**
   * @brief 在指定工作点处构造动力学的一阶线性近似。
   * @param [in] t 当前时间。
   * @param [in] x 线性化状态。
   * @param [in] u 线性化输入。
   * @return 线性近似 `(f, dfdx, dfdu)`。
   */
  VectorFunctionLinearApproximation_t linearApproximation(Scalar t, const StateVector_t& x, const InputVector_t& u) override
  {
    VectorFunctionLinearApproximation_t linearDynamics;
    linearDynamics.f = controlledSystemPtr_->computeFlowMap(t, x, u);
    linearDynamics.dfdx = finiteDifferenceDerivativeState(*controlledSystemPtr_, t, x, u, eps_, doubleSidedDerivative_, isSecondOrderSystem_);
    linearDynamics.dfdu = finiteDifferenceDerivativeInput(*controlledSystemPtr_, t, x, u, eps_, doubleSidedDerivative_, isSecondOrderSystem_);
    return linearDynamics;
  }

private:
  /** @brief 禁用拷贝，避免悬空系统指针被无意复制。 */
  SystemDynamicsLinearizer(const SystemDynamicsLinearizer& other);

  /** @brief 被线性化的非线性系统指针，不拥有对象。 */
  ControlledSystemBase<Scalar, XDimisions, UDimisions>* controlledSystemPtr_;
  /** @brief 为 true 时使用中心差分，否则使用前向差分。 */
  bool doubleSidedDerivative_;
  /** @brief 是否按二阶系统 `[q, q_dot]` 的结构修正 Jacobian。 */
  bool isSecondOrderSystem_;
  /** @brief 基础有限差分步长。 */
  Scalar eps_;
};

