/**
 * @file Types.hpp
 * @brief 基础类型定义：Eigen 矩阵别名与对角矩阵类。
 */
#pragma once
#include <Eigen/Core>
#include <utility>

#include "Matrix/Matrix.hpp"
#include "Matrix/Vector.hpp"

/** @brief 固定尺寸矩阵类型别名，基于 matrix::Matrix。 */
template <typename Scalar, int Rows, int Cols>
using Matrix = matrix::Matrix<Scalar, Rows, Cols>;

/** @brief 固定尺寸列向量类型别名，基于 Eigen::Matrix<Scalar, Rows, 1>。 */
template <typename Scalar, int Rows>
using Vector = matrix::Vector<Scalar, Rows>;

/**
 * @brief 固定维数对角矩阵，仅存储对角线元素，支持与方阵/向量的运算。
 * @tparam Scalar 标量类型（如 double）。
 * @tparam Dimisions 维度（对角线长度）。
 */
template <typename Scalar, int Dimisions>
class DiagonalMatrix {
 public:
  /** @brief 默认构造，元素未初始化。 */
  DiagonalMatrix() = default;
  /** @brief 拷贝构造。 */
  DiagonalMatrix(const DiagonalMatrix& rhs) : data_(rhs.data_) {};
  /** @brief 移动构造。 */
  DiagonalMatrix(DiagonalMatrix&& rhs) : data_(std::move(rhs.data_)) {};

  /**
   * @brief 转换为同尺寸稠密方阵，非对角元素为 0。
   * @return 对角线为 *this 的方阵。
   */
  operator Matrix<Scalar, Dimisions, Dimisions>() {
    Matrix<Scalar, Dimisions, Dimisions> res;
    DiagonalMatrix<Scalar, Dimisions>& self = *this;
    res.setZero();
    for (int i = 0; i < Dimisions; ++i) {
      res(i, i) = self(i);
    }
    return res;
  }

  /** @brief 对角矩阵加法，逐元素。 */
  DiagonalMatrix operator+(const DiagonalMatrix& rhs) const {
    DiagonalMatrix res;
    const DiagonalMatrix& self = *this;
    for (int i = 0; i < Dimisions; ++i) {
      res(i) = self(i) + rhs(i);
    }
    return res;
  }
  /** @brief 对角矩阵复合加法。 */
  DiagonalMatrix& operator+=(const DiagonalMatrix& rhs) {
    DiagonalMatrix& self = *this;
    for (int i = 0; i < Dimisions; ++i) {
      self(i) += rhs(i);
    }
    return self;
  }
  /** @brief 对角矩阵取负。 */
  DiagonalMatrix operator-() const {
    DiagonalMatrix res;
    const DiagonalMatrix& self = *this;
    for (int i = 0; i < Dimisions; ++i) {
      res(i) = -self(i);
    }
    return res;
  }
  /** @brief 对角矩阵减法。 */
  DiagonalMatrix operator-(const DiagonalMatrix& rhs) const {
    DiagonalMatrix res;
    const DiagonalMatrix& self = *this;
    for (int i = 0; i < Dimisions; ++i) {
      res(i) = self(i) - rhs(i);
    }
    return res;
  }
  /** @brief 对角矩阵复合减法。 */
  DiagonalMatrix& operator-=(const DiagonalMatrix& rhs) {
    DiagonalMatrix& self = *this;
    for (int i = 0; i < Dimisions; ++i) {
      self(i) -= rhs(i);
    }
    return self;
  }
  /** @brief 拷贝赋值，含自赋值检查。 */
  DiagonalMatrix& operator=(const DiagonalMatrix& rhs) {
    // 自赋值检查。
    if (this != &rhs) {
      data_ = rhs.data_;
    }
    return *this;
  }
  /** @brief 对角矩阵逐元素乘法。 */
  DiagonalMatrix operator*(const DiagonalMatrix& rhs) const {
    DiagonalMatrix res;
    for (int i = 0; i < Dimisions; ++i) {
      res.data_(i) = data_(i) * rhs.data_(i);
    }
    return res;
  }
  /** @brief 对角矩阵复合逐元素乘法。 */
  DiagonalMatrix& operator*=(const DiagonalMatrix& rhs) {
    for (int i = 0; i < Dimisions; ++i) {
      data_(i) = data_(i) * rhs.data_(i);
    }
    return *this;
  }

  /** @brief 只读访问第 index 个对角元，0 <= index < Dimisions。 */
  Scalar operator()(const int index) const {
    assert(index >= 0 && index < Dimisions);
    return data_(index);
  }
  /** @brief 可写访问第 index 个对角元。 */
  Scalar& operator()(const int index) {
    assert(index >= 0 && index < Dimisions);
    return data_(index);
  }

 private:
  Vector<Scalar, Dimisions> data_;
};

/** @brief 对角矩阵左乘方阵：lhs * rhs，即 rhs 每行乘以 lhs 对应对角元。 */
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator*(
    const DiagonalMatrix<Scalar, Dimisions>& lhs,
    const Matrix<Scalar, Dimisions, Dimisions>& rhs) {
  Matrix<Scalar, Dimisions, Dimisions> res;
  for (int i = 0; i < Dimisions; ++i) {
    res.row(i) = rhs.row(i) * lhs(i);
  }
  return res;
}

/** @brief 方阵右乘对角矩阵：lhs * rhs，即 lhs 每列乘以 rhs 对应对角元。 */
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator*(
    const Matrix<Scalar, Dimisions, Dimisions>& lhs,
    const DiagonalMatrix<Scalar, Dimisions>& rhs) {
  Matrix<Scalar, Dimisions, Dimisions> res;
  for (int i = 0; i < Dimisions; ++i) {
    res.col(i) = lhs.col(i) * rhs(i);
  }
  return res;
}

/** @brief 行向量左乘对角矩阵，逐元素乘。 */
template <typename Scalar, int Dimisions>
Vector<Scalar, Dimisions> operator*(
    const Matrix<Scalar, 1, Dimisions>& lhs,
    const DiagonalMatrix<Scalar, Dimisions>& rhs) {
  Vector<Scalar, Dimisions> res;
  for (int i = 0; i < Dimisions; ++i) {
    res(i) = lhs(i) * rhs(i);
  }
  return res;
}

/** @brief 对角矩阵加方阵：对角元加到方阵对角线上。 */
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator+(
    const DiagonalMatrix<Scalar, Dimisions>& lhs,
    const Matrix<Scalar, Dimisions, Dimisions>& rhs) {
  Matrix<Scalar, Dimisions, Dimisions> res = rhs;
  for (int i = 0; i < Dimisions; ++i) {
    res(i, i) += lhs(i);
  }
  return res;
}
/** @brief 方阵加对角矩阵，委托到 operator+(Diagonal, Matrix)。 */
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator+(
    const Matrix<Scalar, Dimisions, Dimisions>& lhs,
    const DiagonalMatrix<Scalar, Dimisions>& rhs) {
  return rhs + lhs;
}

/** @brief 方阵减对角矩阵：从对角线减去 rhs 对角元。 */
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator-(
    const Matrix<Scalar, Dimisions, Dimisions>& lhs,
    const DiagonalMatrix<Scalar, Dimisions>& rhs) {
  Matrix<Scalar, Dimisions, Dimisions> res = lhs;
  for (int i = 0; i < Dimisions; ++i) {
    res(i, i) -= rhs(i);
  }
  return res;
}

/** @brief 对角矩阵减方阵：结果为 -rhs 再在对角线上加 lhs。 */
template <typename Scalar, int Dimisions>
Matrix<Scalar, Dimisions, Dimisions> operator-(
    const DiagonalMatrix<Scalar, Dimisions>& lhs,
    const Matrix<Scalar, Dimisions, Dimisions>& rhs) {
  Matrix<Scalar, Dimisions, Dimisions> res = -rhs;
  for (int i = 0; i < Dimisions; ++i) {
    res(i, i) += lhs(i);
  }
  return res;
}

template <typename Descriptor>
struct iLQRTypes;