/**
 * @file DiscreteSystemBase.hpp
 * @brief 离散系统基类：仅定义采样节点上的状态转移映射。
 */
#pragma once

#include "Types.hpp"

/**
 * @brief 离散系统基类：子类实现
 * \f$ x_{k+1} = f_d(t_k, x_k, u_k, \Delta t) \f$。
 *
 * 该基类不约定线性化语义；需要雅可比或偏差动力学近似时，应使用
 * DiscreteSystemDynamicsBase。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class DiscreteSystemBase {
 public:
  /** @brief 默认构造。 */
  DiscreteSystemBase() = default;

  /** @brief 虚析构。 */
  virtual ~DiscreteSystemBase() = default;

  /**
   * @brief 给定当前节点的时间、状态、输入和采样周期，计算下一节点状态。
   * @param [in] t 当前采样时刻。
   * @param [in] x 当前状态 \f$ x_k \f$。
   * @param [in] u 当前输入 \f$ u_k \f$。
   * @param [in] dt 当前步长 \f$ t_{k+1} - t_k \f$。
   * @return 下一状态 \f$ x_{k+1} \f$。
   */
  virtual Vector<Scalar, XDim> computeMap(Scalar t,
                                          const Vector<Scalar, XDim>& x,
                                          const Vector<Scalar, UDim>& u,
                                          Scalar dt) = 0;

 protected:
  /** @brief 拷贝构造（供子类使用）。 */
  DiscreteSystemBase(const DiscreteSystemBase& rhs) = default;
};
