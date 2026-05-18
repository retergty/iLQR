/**
 * @file StateConstraint.hpp
 * @brief 仅状态约束基类：约束值及线性/二次近似接口。
 */
#pragma once

#include "ConstraintOrder.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "Types.hpp"
#include <cassert>
#include <cstdlib>

/**
 * @brief
 * 仅状态约束函数基类：按时间与状态返回约束值及线性/二次近似（子类实现）。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 */
template <typename Scalar, int XDim> class StateConstraint {
public:
  /** @brief 构造，指定约束阶数（线性或二次）。 */
  explicit StateConstraint(ConstraintOrder order) : order_(order) {}
  virtual ~StateConstraint() = default;

  /** @brief 获取约束阶数（Linear 或 Quadratic）。 */
  constexpr ConstraintOrder getOrder() const { return order_; };

  /** @brief 获取约束值（标量）。 */
  virtual Scalar getValue(const Scalar time,
                          const Vector<Scalar, XDim> &state) const = 0;

  /** @brief 获取约束的线性近似（默认无效，子类可重写）。 */
  virtual ScalarFunctionLinearApproximation<Scalar, XDim, 0>
  getLinearApproximation(const Scalar time,
                         const Vector<Scalar, XDim> &state) const {
    (void)time;
    (void)state;
    assert(false &&
           "Linear approximation is not implemented for this constraint.");
    std::abort();
  }

  /** @brief 获取约束的二次近似（默认无效，子类可重写）。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>
  getQuadraticApproximation(const Scalar time,
                            const Vector<Scalar, XDim> &state) const {
    (void)time;
    (void)state;
    assert(false &&
           "Quadratic approximation is not implemented for this constraint.");
    std::abort();
  }

private:
  /** @brief 约束阶数。 */
  ConstraintOrder order_;
};
