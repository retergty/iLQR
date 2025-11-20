#pragma once
#include "Dynamics.hpp"

/**
 *
 * A linear time invariant system with the following flow and jump maps:
 *
 * - \f$ \dot{x} = A * x + B * u   \quad \text{for intermediate times}, \f$
 * - \f$ x^{+} = G * x^{-}         \quad \text{for switching times}. \f$
 *
 * where \f$ g(x) \f$ is the guard surface defined by OdeBase::computeGuardSurfaces(t, x).
 */
template <typename Scalar, int XDimisions, int UDimisions>
class LinearSystemDynamics : public SystemDynamicsBase<Scalar, XDimisions, UDimisions>
{
public:
  LinearSystemDynamics(const Matrix<Scalar, XDimisions, XDimisions> &A, const Matrix<Scalar, XDimisions, UDimisions> &B) : A_(A), B_(B) {};

  ~LinearSystemDynamics() override = default;

  Vector<Scalar, XDimisions> computeFlowMap(Scalar t, const Vector<Scalar, XDimisions> &x, const Vector<Scalar, UDimisions> &u) override
  {
    Vector<Scalar, XDimisions> f = A_ * x;
    f += B_ * u;
    return f;
  }

  VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>
  linearApproximation(Scalar t, const Vector<Scalar, XDimisions> &x, const Vector<Scalar, UDimisions> &u) override
  {
    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> approximation;
    approximation.f = A_ * x;
    approximation.f += B_ * u;
    approximation.dfdx = A_;
    approximation.dfdu = B_;
    return approximation;
  }

protected:
  LinearSystemDynamics(const LinearSystemDynamics &other) = default;

  Matrix<Scalar, XDimisions, XDimisions> A_;
  Matrix<Scalar, XDimisions, UDimisions> B_;
};