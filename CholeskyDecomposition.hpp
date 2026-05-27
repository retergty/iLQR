#include "Types.hpp"
#pragma once
// Options 1: solving continuing Problem
template <typename Scalar, int DIMISIONS>
class CholeskyDecomposition {
 public:
  CholeskyDecomposition() = default;
  ~CholeskyDecomposition() = default;
  bool Decomposition(const Matrix<Scalar, DIMISIONS, DIMISIONS>& A) {
    mat_ = A;
    constexpr int n = DIMISIONS;

    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        Scalar sum = mat_(i, j);
        for (int k = i - 1; k >= 0; k--) sum -= mat_(i, k) * mat_(j, k);
        if (i == j) {
          if (sum <= 0) {
            return false;
          }
          mat_(i, i) = std::sqrt(sum);
        } else
          mat_(j, i) = sum / mat_(i, i);
      }
    }
    return true;
  }
  /* Solve L * y = b */
  inline void ForwardElimination(Vector<Scalar, DIMISIONS>& y,
                                 const Vector<Scalar, DIMISIONS>& b) const {
    constexpr int n = DIMISIONS;

    y(0) = b(0) / mat_(0, 0);
    for (int i = 1; i < n; i++) {
      y(i) = b(i);
      for (int j = 0; j < i; j++) y(i) -= mat_(i, j) * y(j);
      y(i) = y(i) / mat_(i, i);
    }
  }
  /* Solve L * Y = B */
  template <int RHS_DIMISIONS>
  inline void ForwardElimination(
      Matrix<Scalar, DIMISIONS, RHS_DIMISIONS>& Y,
      const Matrix<Scalar, DIMISIONS, RHS_DIMISIONS>& B) const {
    constexpr int n = DIMISIONS;

    for (int col = 0; col < RHS_DIMISIONS; ++col) {
      Y(0, col) = B(0, col) / mat_(0, 0);
    }

    for (int i = 1; i < n; i++) {
      for (int col = 0; col < RHS_DIMISIONS; ++col) {
        Y(i, col) = B(i, col);
      }
      for (int j = 0; j < i; j++) {
        for (int col = 0; col < RHS_DIMISIONS; ++col)
          Y(i, col) -= mat_(i, j) * Y(j, col);
      }
      for (int col = 0; col < RHS_DIMISIONS; ++col) {
        Y(i, col) = Y(i, col) / mat_(i, i);
      }
    }
  }
  /* Solve L^T * x = y */
  inline void BackwardElimination(Vector<Scalar, DIMISIONS>& x,
                                  const Vector<Scalar, DIMISIONS>& y) const {
    constexpr int n = DIMISIONS;

    x(n - 1) = y(n - 1) / mat_(n - 1, n - 1);
    for (int i = n - 2; i >= 0; i--) {
      x(i) = y(i);
      for (int j = i + 1; j < n; j++) x(i) -= mat_(j, i) * x(j);
      x(i) = x(i) / mat_(i, i);
    }
  }
  /* Solve L^T * X = Y */
  template <int RHS_DIMISIONS>
  inline void BackwardElimination(
      Matrix<Scalar, DIMISIONS, RHS_DIMISIONS>& X,
      const Matrix<Scalar, DIMISIONS, RHS_DIMISIONS>& Y) const {
    constexpr int n = DIMISIONS;

    for (int col = 0; col < RHS_DIMISIONS; ++col) {
      X(n - 1, col) = Y(n - 1, col) / mat_(n - 1, n - 1);
    }

    for (int i = n - 2; i >= 0; i--) {
      for (int col = 0; col < RHS_DIMISIONS; ++col) {
        X(i, col) = Y(i, col);
      }
      for (int j = i + 1; j < n; j++) {
        for (int col = 0; col < RHS_DIMISIONS; ++col)
          X(i, col) -= mat_(j, i) * X(j, col);
      }
      for (int col = 0; col < RHS_DIMISIONS; ++col) {
        X(i, col) = X(i, col) / mat_(i, i);
      }
    }
  }
  // solve L*L^T*x = b
  void Solve(Vector<Scalar, DIMISIONS>& x,
             const Vector<Scalar, DIMISIONS>& b) const {
    Vector<Scalar, DIMISIONS> y;

    /* Solve L * y = b */
    ForwardElimination(y, b);
    /* Solve L^T * x = y */
    BackwardElimination(x, y);
  }
  // solve L*L^T*X = B
  template <int RHS_DIMISIONS>
  void Solve(Matrix<Scalar, DIMISIONS, RHS_DIMISIONS>& X,
             const Matrix<Scalar, DIMISIONS, RHS_DIMISIONS>& B) const {
    Matrix<Scalar, DIMISIONS, RHS_DIMISIONS> Y;

    /* Solve L * Y = B */
    ForwardElimination(Y, B);
    /* Solve L^T * X = Y */
    BackwardElimination(X, Y);
  }

  Matrix<Scalar, DIMISIONS, DIMISIONS> GetDecompositionResult() const {
    Matrix<Scalar, DIMISIONS, DIMISIONS> result;
    for (int i = 0; i < DIMISIONS; ++i) {
      result(i, i) = mat_(i, i);
      for (int j = 0; j < i; ++j) {
        result(i, j) = mat_(i, j);
        result(j, i) = mat_(i, j);
      }
    }
    return result;
  }

  Matrix<Scalar, DIMISIONS, DIMISIONS> GetMatrixL() const {
    Matrix<Scalar, DIMISIONS, DIMISIONS> L;

    for (int i = 0; i < DIMISIONS; ++i) {
      for (int j = 0; j <= i; ++j) {
        L(i, j) = mat_(i, j);
      }
      for (int j = i + 1; j < DIMISIONS; ++j) {
        L(i, j) = 0;
      }
    }
    return L;
  }
  Matrix<Scalar, DIMISIONS, DIMISIONS> GetMatrixLT() const {
    Matrix<Scalar, DIMISIONS, DIMISIONS> LT;

    for (int i = 0; i < DIMISIONS; ++i) {
      for (int j = 0; j < i; ++j) {
        LT(i, j) = 0;
      }
      for (int j = i; j < DIMISIONS; ++j) {
        LT(i, j) = mat_(j, i);
      }
    }
    return LT;
  }
  Matrix<Scalar, DIMISIONS, DIMISIONS> InverseL() const {
    Matrix<Scalar, DIMISIONS, DIMISIONS> inv_L;
    Vector<Scalar, DIMISIONS> y;
    Vector<Scalar, DIMISIONS> b;
    b.setZero();
    // L^-1
    for (int i = 0; i < DIMISIONS; ++i) {
      b(i) = 1;
      this->ForwardElimination(y, b);
      inv_L.col(i) = y;
      b(i) = 0;
    }
    return inv_L;
  }
  Matrix<Scalar, DIMISIONS, DIMISIONS> Inverse() const {
    Matrix<Scalar, DIMISIONS, DIMISIONS> inv_L = InverseL();

    Matrix<Scalar, DIMISIONS, DIMISIONS> inv;
    Vector<Scalar, DIMISIONS> x;
    Vector<Scalar, DIMISIONS> y;
    // L^T*X=Y
    for (int i = 0; i < DIMISIONS; ++i) {
      y = inv_L.col(i);
      this->BackwardElimination(x, y);
      inv.col(i) = x;
    }

    return inv;
  }

 private:
  Matrix<Scalar, DIMISIONS, DIMISIONS> mat_;  // L*L^T decompsition result
};