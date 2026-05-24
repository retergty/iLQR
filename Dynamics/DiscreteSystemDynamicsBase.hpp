/**
 * @file DiscreteSystemDynamicsBase.hpp
 * @brief 离散系统动力学与线性化基类。
 *
 * 离散状态转移：\f$ x_{k+1} = f_d(t_k, x_k, u_k, \Delta t) \f$
 *
 * 完整一阶近似：
 * \f[
 *   x_{k+1} \approx x_{next,nom} + A_k\delta x_k + B_k\delta u_k
 * \f]
 *
 * 偏差动力学近似：
 * \f[
 *   \delta x_{k+1} \approx A_k\delta x_k + B_k\delta u_k
 * \f]
 */
#pragma once

#include "DiscreteSystemBase.hpp"
#include "LinearApproximation.hpp"

/**
 * @brief 离散系统动力学基类：子类需实现离散状态转移及其线性近似。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class DiscreteSystemDynamicsBase
    : public DiscreteSystemBase<Scalar, XDim, UDim> {
 public:
  using LinearApproximation_t =
      VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>;

  /** @brief 默认构造。 */
  DiscreteSystemDynamicsBase() = default;

  /** @brief 虚析构。 */
  ~DiscreteSystemDynamicsBase() override = default;

  /**
   * @brief 计算离散状态转移的完整一阶近似。
   * @param [in] t 当前采样时刻。
   * @param [in] x 当前状态 \f$ x_k \f$。
   * @param [in] u 当前输入 \f$ u_k \f$。
   * @param [in] dt 当前步长 \f$ t_{k+1} - t_k \f$。
   * @return 包含名义下一状态、dfdx 和 dfdu 的离散一阶近似。
   *
   * 返回值中的 dfdx 和 dfdu 表示离散映射
   * \f$ f_d(t_k,x_k,u_k,\Delta t) \f$ 对状态和输入的雅可比。
   * f 表示名义下一状态。
   */
  virtual LinearApproximation_t linearApproximation(
      Scalar t, const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u,
      Scalar dt) = 0;

  /**
   * @brief 计算离散偏差动力学近似，仅保留 dfdx、dfdu 并将常数项置零。
   *
   * 默认实现复用完整线性化以保持旧模型兼容。性能敏感模型可重写该函数，
   * 直接填充雅可比并避免计算名义下一状态。
   */
  virtual LinearApproximation_t deviationLinearApproximation(
      Scalar t, const Vector<Scalar, XDim>& x, const Vector<Scalar, UDim>& u,
      Scalar dt) {
    LinearApproximation_t approximation = linearApproximation(t, x, u, dt);
    approximation.f.setZero();
    return approximation;
  }

 protected:
  /** @brief 拷贝构造（供子类使用）。 */
  DiscreteSystemDynamicsBase(const DiscreteSystemDynamicsBase& rhs) = default;
};
