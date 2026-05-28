/**
 * @file SystemDynamicsBase.hpp
 * @brief 连续时间系统动力学与线性化基类。
 *
 * 子类实现连续流映射 \f$ \dot{x}=f(t,x,u) \f$ 及其一阶近似：
 * \f[
 *   \dot{x} \approx f_{nom} + A(t)\delta x + B(t)\delta u
 * \f]
 *
 * 求解器若只需要偏差动力学：
 * \f[
 *   \delta\dot{x} \approx A(t)\delta x + B(t)\delta u
 * \f]
 * 应使用 deviationLinearApproximation()，而不是由本接口隐含 Riccati
 * 前的常数项处理。
 */
#pragma once
#include "Approximation/LinearApproximation.hpp"
#include "Dynamics/ControlledSystemBase.hpp"

/**
 * @brief 系统动力学与线性化基类：子类需实现 linearApproximation。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class SystemDynamicsBase : public ControlledSystemBase<Scalar, XDim, UDim> {
 public:
  using LinearApproximation_t =
      VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>;

  /** @brief 默认构造。 */
  SystemDynamicsBase() = default;

  /** @brief 虚析构函数。 */
  ~SystemDynamicsBase() override = default;

  /**
   * @brief 计算工作点处连续动力学的完整一阶近似。
   * @param [in] t 当前时间。
   * @param [in] x 当前状态。
   * @param [in] u 当前输入。
   * @return 包含 \f$f_{nom}=f(t,x,u)\f$、dfdx 和 dfdu 的线性近似。
   */
  virtual LinearApproximation_t linearApproximation(
      Scalar t, const Vector<Scalar, XDim>& x,
      const Vector<Scalar, UDim>& u) = 0;

  /**
   * @brief 计算偏差动力学近似，仅保留 dfdx、dfdu 并将常数项置零。
   *
   * 默认实现复用完整线性化以保持旧模型兼容。性能敏感模型可重写该函数，
   * 直接填充雅可比并避免计算 \f$f_{nom}\f$。
   */
  virtual LinearApproximation_t deviationLinearApproximation(
      Scalar t, const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u) {
    LinearApproximation_t approximation = linearApproximation(t, x, u);
    approximation.f.setZero();
    return approximation;
  }
};
