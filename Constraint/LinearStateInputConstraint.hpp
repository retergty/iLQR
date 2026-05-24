#pragma once

#include "StateInputConstraint.hpp"
#include "Types.hpp"

/**
 * 线性状态-输入约束。
 */
template <typename Scalar, int XDim, int UDim, int CDim>
class LinearStateInputConstraint
    : public StateInputConstraint<Scalar, XDim, UDim, CDim> {
 public:
  /**
   * 构造函数
   *
   * @param[in] e: C * x + D * u + e = 0 中的常数项。
   * @param[in] C: C * x + D * u + e = 0 中的 x 系数。
   * @param[in] D: C * x + D * u + e = 0 中的 u 系数。
   */
  LinearStateInputConstraint(const Vector<Scalar, CDim>& e,
                             const Matrix<Scalar, CDim, XDim>& C,
                             const Matrix<Scalar, CDim, UDim>& D)
      : StateInputConstraint<Scalar, XDim, UDim, CDim>(ConstraintOrder::Linear),
        e_(e),
        C_(C),
        D_(D) {}

  ~LinearStateInputConstraint() override = default;

  Vector<Scalar, CDim> getValue(const Scalar time,
                                const Vector<Scalar, XDim>& state,
                                const Vector<Scalar, UDim>& input) const final {
    (void)time;
    Vector<Scalar, CDim> g = e_;
    g += C_ * state;
    g += D_ * input;
    return g;
  }

  VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim>
  getLinearApproximation(const Scalar time, const Vector<Scalar, XDim>& state,
                         const Vector<Scalar, UDim>& input) const final {
    (void)time;
    VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim> g;
    g.f = e_;
    g.f += C_ * state;
    g.f += D_ * input;
    g.dfdx = C_;
    g.dfdu = D_;
    return g;
  }

 public:
  Vector<Scalar, CDim> e_;       /**< 状态-输入约束。 */
  Matrix<Scalar, CDim, XDim> C_; /**< 状态-输入约束对状态的导数。 */
  Matrix<Scalar, CDim, UDim> D_; /**< 状态-输入约束对输入的导数。 */
};
