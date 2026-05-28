#pragma once
#include "Constraint/StateConstraint.hpp"
#include "Types.hpp"

/**
 * 线性仅状态约束。
 */
template <typename Scalar, int XDim, int CDim>
class LinearStateConstraint : public StateConstraint<Scalar, XDim, CDim> {
 public:
  /**
   * 构造函数
   *
   * @param [in] h: F * x + h = 0 中的常数项。
   * @param [in] F: F * x + h = 0 中的 x 系数。
   */
  LinearStateConstraint(const Vector<Scalar, CDim>& h,
                        const Matrix<Scalar, CDim, XDim>& F)
      : StateConstraint<Scalar, XDim, CDim>(ConstraintOrder::Linear),
        h_(h),
        F_(F) {}

  ~LinearStateConstraint() override = default;

  Vector<Scalar, CDim> getValue(const Scalar time,
                                const Vector<Scalar, XDim>& state) const final {
    (void)time;
    Vector<Scalar, CDim> g = h_;
    g += F_ * state;
    return g;
  }

  VectorFunctionLinearApproximation<Scalar, CDim, XDim, 0>
  getLinearApproximation(const Scalar time,
                         const Vector<Scalar, XDim>& state) const final {
    (void)time;
    VectorFunctionLinearApproximation<Scalar, CDim, XDim, 0> g;
    g.f = h_;
    g.f += F_ * state;
    g.dfdx = F_;
    return g;
  }

 public:
  Vector<Scalar, CDim> h_;       /**< 仅状态约束。 */
  Matrix<Scalar, CDim, XDim> F_; /**< 仅状态约束对状态的导数。 */
};
