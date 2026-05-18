/**
 * @file StateInputConstraint.hpp
 * @brief 状态-输入约束基类：约束值及线性/二次近似接口。
 */
#pragma once
#include <cassert>
#include <cstdlib>

#include "ConstraintOrder.hpp"
#include "IntrusiveList.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "Types.hpp"

/**
 * @brief
 * 状态-输入约束函数基类：按时间、状态与输入返回约束值及线性/二次近似（子类实现）。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int XDim, int UDim>
class StateInputConstraint {
 public:
  /** @brief 构造，指定约束阶数（线性或二次）。 */
  explicit StateInputConstraint(ConstraintOrder order) : order_(order) {}
  virtual ~StateInputConstraint() = default;

  /** @brief 获取约束阶数（Linear 或 Quadratic）。 */
  constexpr ConstraintOrder getOrder() const { return order_; };

  /** @brief 获取约束值（标量）。 */
  virtual Scalar getValue(const Scalar time, const Vector<Scalar, XDim>& state,
                          const Vector<Scalar, UDim>& input) const = 0;

  /** @brief 获取约束的线性近似（默认无效，子类可重写）。 */
  virtual ScalarFunctionLinearApproximation<Scalar, XDim, UDim>
  getLinearApproximation(const Scalar time, const Vector<Scalar, XDim>& state,
                         const Vector<Scalar, UDim>& input) const {
    (void)time;
    (void)state;
    (void)input;
    assert(false &&
           "Linear approximation is not implemented for this constraint.");
    return ScalarFunctionLinearApproximation<Scalar, XDim, 0>::Zero();
  }

  /** @brief 获取约束的二次近似（默认无效，子类可重写）。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(const Scalar time,
                            const Vector<Scalar, XDim>& state,
                            const Vector<Scalar, UDim>& input) const {
    (void)time;
    (void)state;
    (void)input;
    assert(false &&
           "Quadratic approximation is not implemented for this constraint.");
    return ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>::Zero();
  }

 private:
  /** @brief 约束阶数。 */
  ConstraintOrder order_;
};