/**
 * @file matrix/Types.hpp
 * @brief 基础固定尺寸矩阵与向量类型别名。
 */
#pragma once

#include "matrix/Matrix.hpp"
#include "matrix/Vector.hpp"

/** @brief 固定尺寸矩阵类型别名，基于 matrix::Matrix。 */
template <typename Scalar, int Rows, int Cols>
using Matrix = matrix::Matrix<Scalar, Rows, Cols>;

/** @brief 固定尺寸列向量类型别名，基于 matrix::Vector。 */
template <typename Scalar, int Rows>
using Vector = matrix::Vector<Scalar, Rows>;
