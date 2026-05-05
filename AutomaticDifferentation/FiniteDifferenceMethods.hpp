#pragma once
#include "Types.hpp"
#include <functional>
#include <math.h>
#include "ControlledSystemBase.hpp"

template<typename Scalar, int StateDimisions, int VarDimisions>
Matrix<Scalar, StateDimisions, VarDimisions> finiteDifferenceDerivative(std::function<Vector<Scalar, StateDimisions>(const Vector<Scalar, VarDimisions>&)> f,
  const Vector<Scalar, VarDimisions>& x0, Scalar eps,
  bool doubleSidedDerivative) {
  const Vector<Scalar, StateDimisions> f0 = f(x0);
  Matrix<Scalar, StateDimisions, VarDimisions> jacobian;

  for (size_t i = 0; i < VarDimisions; i++) {
    // inspired from: http://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating_point_arithmetic
    Scalar h = eps * std::max(std::fabs(x0(i)), 1.0);

    Vector<Scalar, VarDimisions> xPlus = x0;
    xPlus(i) += h;

    if (doubleSidedDerivative) {
      Vector<Scalar, VarDimisions> xMinus = x0;
      xMinus(i) -= h;
      jacobian.col(i) = (f(xPlus) - f(xMinus)) / (2.0 * h);
    }
    else {
      jacobian.col(i) = (f(xPlus) - f0) / h;
    }
  }

  return jacobian;
}


template<typename Scalar, int XDimisions, int UDimisions>
Matrix<Scalar, XDimisions, XDimisions> finiteDifferenceDerivativeState(ControlledSystemBase<Scalar, XDimisions, UDimisions>& system,
  Scalar t, const Vector<Scalar, XDimisions>& x, const  Vector<Scalar, UDimisions>& u, Scalar eps,
  bool doubleSidedDerivative, bool isSecondOrderSystem) {
  auto f = [&](const Vector<Scalar, XDimisions>& var) -> Vector<Scalar, XDimisions> { return system.computeFlowMap(t, var, u); };

  Matrix<Scalar, XDimisions, XDimisions> A = finiteDifferenceDerivative(f, x, eps, doubleSidedDerivative);

  if (isSecondOrderSystem) {
    // Assumes state vector = [x, x_dot]
    A.topLeftCorner(x.rows() / 2, x.rows() / 2).setZero();
    A.topRightCorner(x.rows() / 2, x.rows() / 2).setIdentity();
  }
  return A;
}

template<typename Scalar, int XDimisions, int UDimisions>
Matrix<Scalar, XDimisions, UDimisions> finiteDifferenceDerivativeInput(ControlledSystemBase<Scalar, XDimisions, UDimisions>& system,
  Scalar t, const Vector<Scalar, XDimisions>& x, const  Vector<Scalar, UDimisions>& u, Scalar eps,
  bool doubleSidedDerivative, bool isSecondOrderSystem)
{
  auto f = [&](const Vector<Scalar, UDimisions>& var) -> Vector<Scalar, XDimisions> { return system.computeFlowMap(t, x, var); };

  Matrix<Scalar, XDimisions, UDimisions> B = finiteDifferenceDerivative(f, u, eps, doubleSidedDerivative);

  if (isSecondOrderSystem) {
    // Assumes state vector = [x, x_dot]
    B.topRows(x.rows() / 2).setZero();
  }
  return B;
}