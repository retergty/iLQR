#pragma once

#include <Eigen/Core>
#include <cassert>

#include "matrix/Matrix.hpp"
#include "matrix/Vector.hpp"

namespace test_tools {
namespace matrix_eigen_conversion {
namespace detail {

template <int Dim>
struct EigenDim {
  static_assert(Dim >= 0, "Custom matrix dimension must be non-negative.");
  static constexpr int value = Dim;
};

template <int Dim>
struct FixedEigenDim {
  static_assert(Dim != Eigen::Dynamic,
                "Only fixed-size Eigen objects can convert to matrix::Matrix.");
  static_assert(Dim >= 0, "Eigen dimension must be non-negative.");
  static constexpr int value = Dim;
};

}  // namespace detail

template <typename Scalar, int Rows, int Cols>
using EigenMatrix_t = Eigen::Matrix<Scalar, detail::EigenDim<Rows>::value,
                                    detail::EigenDim<Cols>::value>;

template <typename Scalar, int Rows>
using EigenVector_t = Eigen::Matrix<Scalar, detail::EigenDim<Rows>::value, 1>;

template <typename Scalar, int Rows, int Cols>
EigenMatrix_t<Scalar, Rows, Cols> toEigenMatrix(
    const matrix::Matrix<Scalar, Rows, Cols>& src) {
  EigenMatrix_t<Scalar, Rows, Cols> dst;
  for (int i = 0; i < Rows; ++i) {
    for (int j = 0; j < Cols; ++j) {
      dst(i, j) = src(i, j);
    }
  }
  return dst;
}

template <typename Scalar, int Rows>
EigenVector_t<Scalar, Rows> toEigenVector(
    const matrix::Vector<Scalar, Rows>& src) {
  EigenVector_t<Scalar, Rows> dst;
  for (int i = 0; i < Rows; ++i) {
    dst(i) = src(i);
  }
  return dst;
}

template <typename Derived>
matrix::Matrix<typename Derived::Scalar,
               detail::FixedEigenDim<Derived::RowsAtCompileTime>::value,
               detail::FixedEigenDim<Derived::ColsAtCompileTime>::value>
fromEigenMatrix(const Eigen::MatrixBase<Derived>& src) {
  using Scalar = typename Derived::Scalar;
  constexpr int CustomRows =
      detail::FixedEigenDim<Derived::RowsAtCompileTime>::value;
  constexpr int CustomCols =
      detail::FixedEigenDim<Derived::ColsAtCompileTime>::value;

  assert(src.rows() == CustomRows);
  assert(src.cols() == CustomCols);

  matrix::Matrix<Scalar, CustomRows, CustomCols> dst;
  for (int i = 0; i < CustomRows; ++i) {
    for (int j = 0; j < CustomCols; ++j) {
      dst(i, j) = src(i, j);
    }
  }
  return dst;
}

template <int Rows, int Cols, typename Derived>
matrix::Matrix<typename Derived::Scalar, Rows, Cols> fromEigenMatrix(
    const Eigen::MatrixBase<Derived>& src) {
  using Scalar = typename Derived::Scalar;

  assert(src.rows() == Rows);
  assert(src.cols() == Cols);

  matrix::Matrix<Scalar, Rows, Cols> dst;
  for (int i = 0; i < Rows; ++i) {
    for (int j = 0; j < Cols; ++j) {
      dst(i, j) = src(i, j);
    }
  }
  return dst;
}

template <typename Derived>
matrix::Vector<typename Derived::Scalar,
               detail::FixedEigenDim<Derived::RowsAtCompileTime>::value>
fromEigenVector(const Eigen::MatrixBase<Derived>& src) {
  static_assert(Derived::ColsAtCompileTime == 1,
                "fromEigenVector requires a column vector expression.");

  using Scalar = typename Derived::Scalar;
  constexpr int CustomRows =
      detail::FixedEigenDim<Derived::RowsAtCompileTime>::value;

  assert(src.rows() == CustomRows);
  assert(src.cols() == 1);

  matrix::Vector<Scalar, CustomRows> dst;
  for (int i = 0; i < CustomRows; ++i) {
    dst(i) = src(i);
  }
  return dst;
}

template <int Rows, typename Derived>
matrix::Vector<typename Derived::Scalar, Rows> fromEigenVector(
    const Eigen::MatrixBase<Derived>& src) {
  using Scalar = typename Derived::Scalar;

  assert(src.rows() == Rows);
  assert(src.cols() == 1);

  matrix::Vector<Scalar, Rows> dst;
  for (int i = 0; i < Rows; ++i) {
    dst(i) = src(i);
  }
  return dst;
}

}  // namespace matrix_eigen_conversion
}  // namespace test_tools
