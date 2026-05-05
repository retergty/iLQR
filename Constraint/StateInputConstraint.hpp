/**
 * @file StateInputConstraint.hpp
 * @brief 状态-输入约束基类：约束值及线性/二次近似接口。
 */
#pragma once
#include "Types.hpp"
#include "ConstraintOrder.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "IntrusiveList.hpp"

/**
 * @brief 状态-输入约束函数基类：按时间、状态与输入返回约束值及线性/二次近似（子类实现）。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam UDimisions 输入维度。
 */
template<typename Scalar, int XDimisions, int UDimisions>
class StateInputConstraint
{
public:
  /** @brief 构造，指定约束阶数（线性或二次）。 */
  explicit StateInputConstraint(ConstraintOrder order) : order_(order) {}
  virtual ~StateInputConstraint() = default;

  /** @brief 获取约束阶数（Linear 或 Quadratic）。 */
  constexpr ConstraintOrder getOrder() const { return order_; };

  /** @brief 获取约束值（标量）。 */
  virtual Scalar getValue(const Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input) const = 0;

  /** @brief 获取约束的线性近似（默认无效，子类可重写）。 */
  virtual ScalarFunctionLinearApproximation<Scalar, XDimisions, UDimisions> getLinearApproximation(const Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input) const
  {
    static_assert(!std::is_same_v<Scalar, Scalar>, "invalid!");
  }

  /** @brief 获取约束的二次近似（默认无效，子类可重写）。 */
  virtual ScalarFunctionLinearApproximation<Scalar, XDimisions, UDimisions> getQuadraticApproximation(const Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input) const {
    static_assert(!std::is_same_v<Scalar, Scalar>, "invalid!");
  }

private:
  /** @brief 约束阶数。 */
  ConstraintOrder order_;
};