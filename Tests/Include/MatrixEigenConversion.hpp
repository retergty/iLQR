#pragma once

#include <Eigen/Core>

#include <cassert>
#include <cstddef>
#include <limits>

#include "../../Matrix/Matrix.hpp"
#include "../../Matrix/Vector.hpp"

namespace test_tools {
namespace matrix_eigen_conversion {
namespace detail {

template <std::size_t Dim>
struct EigenDim {
  static_assert(Dim <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
                "Custom matrix dimension exceeds Eigen's int dimension range.");
  static constexpr int value = static_cast<int>(Dim);
};

template <int Dim>
struct FixedEigenDim {
  static_assert(Dim != Eigen::Dynamic,
                "Only fixed-size Eigen objects can convert to matrix::Matrix.");
  static_assert(Dim >= 0, "Eigen dimension must be non-negative.");
  static constexpr std::size_t value = static_cast<std::size_t>(Dim);
};

}  // namespace detail

template <typename Scalar, std::size_t Rows, std::size_t Cols>
using EigenMatrix_t =
    Eigen::Matrix<Scalar, detail::EigenDim<Rows>::value,
                  detail::EigenDim<Cols>::value>;

template <typename Scalar, std::size_t Rows>
using EigenVector_t =
    Eigen::Matrix<Scalar, detail::EigenDim<Rows>::value, 1>;

template <typename Scalar, std::size_t Rows, std::size_t Cols>
EigenMatrix_t<Scalar, Rows, Cols> toEigenMatrix(
    const matrix::Matrix<Scalar, Rows, Cols>& src) {
  EigenMatrix_t<Scalar, Rows, Cols> dst;
  for (std::size_t i = 0; i < Rows; ++i) {
    for (std::size_t j = 0; j < Cols; ++j) {
      dst(static_cast<int>(i), static_cast<int>(j)) = src(i, j);
    }
  }
  return dst;
}

template <typename Scalar, std::size_t Rows>
EigenVector_t<Scalar, Rows> toEigenVector(
    const matrix::Vector<Scalar, Rows>& src) {
  EigenVector_t<Scalar, Rows> dst;
  for (std::size_t i = 0; i < Rows; ++i) {
    dst(static_cast<int>(i)) = src(i);
  }
  return dst;
}

template <typename Derived>
matrix::Matrix<typename Derived::Scalar,
               detail::FixedEigenDim<Derived::RowsAtCompileTime>::value,
               detail::FixedEigenDim<Derived::ColsAtCompileTime>::value>
fromEigenMatrix(const Eigen::MatrixBase<Derived>& src) {
  using Scalar = typename Derived::Scalar;
  constexpr std::size_t CustomRows =
      detail::FixedEigenDim<Derived::RowsAtCompileTime>::value;
  constexpr std::size_t CustomCols =
      detail::FixedEigenDim<Derived::ColsAtCompileTime>::value;

  assert(src.rows() == static_cast<int>(CustomRows));
  assert(src.cols() == static_cast<int>(CustomCols));

  matrix::Matrix<Scalar, CustomRows, CustomCols> dst;
  for (std::size_t i = 0; i < CustomRows; ++i) {
    for (std::size_t j = 0; j < CustomCols; ++j) {
      dst(i, j) = src(static_cast<int>(i), static_cast<int>(j));
    }
  }
  return dst;
}

template <std::size_t Rows, std::size_t Cols, typename Derived>
matrix::Matrix<typename Derived::Scalar, Rows, Cols> fromEigenMatrix(
    const Eigen::MatrixBase<Derived>& src) {
  using Scalar = typename Derived::Scalar;

  assert(src.rows() == static_cast<int>(Rows));
  assert(src.cols() == static_cast<int>(Cols));

  matrix::Matrix<Scalar, Rows, Cols> dst;
  for (std::size_t i = 0; i < Rows; ++i) {
    for (std::size_t j = 0; j < Cols; ++j) {
      dst(i, j) = src(static_cast<int>(i), static_cast<int>(j));
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
  constexpr std::size_t CustomRows =
      detail::FixedEigenDim<Derived::RowsAtCompileTime>::value;

  assert(src.rows() == static_cast<int>(CustomRows));
  assert(src.cols() == 1);

  matrix::Vector<Scalar, CustomRows> dst;
  for (std::size_t i = 0; i < CustomRows; ++i) {
    dst(i) = src(static_cast<int>(i));
  }
  return dst;
}

template <std::size_t Rows, typename Derived>
matrix::Vector<typename Derived::Scalar, Rows> fromEigenVector(
    const Eigen::MatrixBase<Derived>& src) {
  using Scalar = typename Derived::Scalar;

  assert(src.rows() == static_cast<int>(Rows));
  assert(src.cols() == 1);

  matrix::Vector<Scalar, Rows> dst;
  for (std::size_t i = 0; i < Rows; ++i) {
    dst(i) = src(static_cast<int>(i));
  }
  return dst;
}

}  // namespace matrix_eigen_conversion
}  // namespace test_tools
