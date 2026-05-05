#pragma once
#include "ControlledSystemBase.hpp"
#include "SystemDynamicsBase.hpp"
/**
 * @file LinearSystemDynamics.hpp
 * @brief 线性时不变系统动力学：流映射 \f$ \dot{x} = A x + B u \f$。
 */
template <typename Scalar, int XDimisions, int UDimisions>
class LinearSystemDynamics : public SystemDynamicsBase<Scalar, XDimisions, UDimisions>
{
public:
    using SystemDynamicsBase<Scalar, XDimisions, UDimisions>::computeFlowMap;
  /**
   * @brief 用状态矩阵 A 与控制矩阵 B 构造线性系统。
   * @param [in] A 状态矩阵。
   * @param [in] B 控制矩阵。
   */
  LinearSystemDynamics(const Matrix<Scalar, XDimisions, XDimisions> &A, const Matrix<Scalar, XDimisions, UDimisions> &B) : A_(A), B_(B) {};

  /** @brief 析构函数。 */
  ~LinearSystemDynamics() override = default;

  /** @brief 计算流映射：\f$ \dot{x} = A x + B u \f$。 */
  Vector<Scalar, XDimisions> computeFlowMap(Scalar t, const Vector<Scalar, XDimisions> &x, const Vector<Scalar, UDimisions> &u) const override
  {
    (void)t;
    Vector<Scalar, XDimisions> f = A_ * x;
    f += B_ * u;
    return f;
  }

  /** @brief 返回线性近似：dfdx=A, dfdu=B, f=A*x+B*u。 */
  VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>
  linearApproximation(Scalar t, const Vector<Scalar, XDimisions> &x, const Vector<Scalar, UDimisions> &u) override
  {
    (void)t;
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> approximation;
    approximation.f = A_ * x;
    approximation.f += B_ * u;
    approximation.dfdx = A_;
    approximation.dfdu = B_;
    return approximation;
  }

protected:
  /** @brief 拷贝构造（供子类使用）。 */
  LinearSystemDynamics(const LinearSystemDynamics &other) = default;

  /** @brief 状态矩阵 A。 */
  Matrix<Scalar, XDimisions, XDimisions> A_;
  /** @brief 控制矩阵 B。 */
  Matrix<Scalar, XDimisions, UDimisions> B_;
};