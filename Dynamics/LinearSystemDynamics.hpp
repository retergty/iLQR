#pragma once
#include "ControlledSystemBase.hpp"
#include "SystemDynamicsBase.hpp"
/**
 * @file LinearSystemDynamics.hpp
 * @brief 线性时不变系统动力学：流映射 \f$ \dot{x} = A x + B u \f$。
 */
template <typename Scalar, int XDim, int UDim>
class LinearSystemDynamics : public SystemDynamicsBase<Scalar, XDim, UDim>
{
public:
    using SystemDynamicsBase<Scalar, XDim, UDim>::computeFlowMap;
  /**
   * @brief 用状态矩阵 A 与控制矩阵 B 构造线性系统。
   * @param [in] A 状态矩阵。
   * @param [in] B 控制矩阵。
   */
  LinearSystemDynamics(const Matrix<Scalar, XDim, XDim> &A, const Matrix<Scalar, XDim, UDim> &B) : A_(A), B_(B) {};

  /** @brief 析构函数。 */
  ~LinearSystemDynamics() override = default;

  /** @brief 计算流映射：\f$ \dot{x} = A x + B u \f$。 */
  Vector<Scalar, XDim> computeFlowMap(Scalar t, const Vector<Scalar, XDim> &x, const Vector<Scalar, UDim> &u) const override
  {
    (void)t;
    Vector<Scalar, XDim> f = A_ * x;
    f += B_ * u;
    return f;
  }

  /** @brief 返回线性近似：dfdx=A, dfdu=B, f=A*x+B*u。 */
  VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
  linearApproximation(Scalar t, const Vector<Scalar, XDim> &x, const Vector<Scalar, UDim> &u) override
  {
    (void)t;
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> approximation;
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
  Matrix<Scalar, XDim, XDim> A_;
  /** @brief 控制矩阵 B。 */
  Matrix<Scalar, XDim, UDim> B_;
};