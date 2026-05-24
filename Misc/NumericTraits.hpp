/**
 * @file NumericTraits.hpp
 * @brief 数值特征：比较与极限判定的 epsilon 常量（limitEpsilon、weakEpsilon）。
 */
#pragma once

namespace numeric_traits {
/**
 * @brief 极限邻域 epsilon：若 v 在 w 的该邻域内则视为 v 趋于 w，应大于
 * weakEpsilon。
 * @return limit epsilon 值（默认 1e-6）。
 */
template <typename T>
constexpr T limitEpsilon() {
  return T(1e-6);
}

/**
 * @brief 比较时使用的弱精度 epsilon。
 * @return weak epsilon 值（默认 1e-9）。
 */
template <typename T>
constexpr T weakEpsilon() {
  return T(1e-9);
}

}  // namespace numeric_traits