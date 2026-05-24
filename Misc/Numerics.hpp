/**
 * @file Numerics.hpp
 * @brief 数值比较工具：基于机器精度的浮点近似相等、近似小于等于、近似大于等于。
 */
#pragma once

#include <cmath>
#include <limits>
#include <type_traits>

namespace numerics {
/**
 * @brief 在给定精度下判断两浮点数是否近似相等（考虑量级与机器精度）。
 * @tparam T1 第一个参数类型。
 * @tparam T2 第二个参数类型。
 * @tparam T3 精度类型。
 * @param [in] x 第一个浮点数。
 * @param [in] y 第二个浮点数。
 * @param [in] prec 比较精度。
 * @return 近似相等返回 true。
 */
template <class T1, class T2, class T3>
constexpr bool almost_eq(T1&& x, T2&& y, T3&& prec) {
  static_assert(
      std::is_floating_point<typename std::remove_reference<T1>::type>::value,
      "First argument is not floating point!");
  static_assert(
      std::is_floating_point<typename std::remove_reference<T2>::type>::value,
      "Second argument is not floating point!");
  static_assert(
      std::is_floating_point<typename std::remove_reference<T3>::type>::value,
      "prec is not floating point!");
  // 机器 epsilon 需要按所用数值的量级缩放，
  // 并乘以期望精度，除非结果为次正规数。
  using Type = const std::remove_reference_t<T1>;
  const auto absDiff = std::abs(x - static_cast<Type>(y));
  const auto magnitude = std::min(std::abs(x), std::abs(static_cast<Type>(y)));
  return absDiff <= static_cast<Type>(prec) * magnitude ||
         absDiff < std::numeric_limits<Type>::min();
}

/**
 * @brief 使用机器精度判断两浮点数是否近似相等。
 * @param [in] x 第一个浮点数。
 * @param [in] y 第二个浮点数。
 * @return 近似相等返回 true。
 */
template <class T1, class T2>
constexpr bool almost_eq(T1&& x, T2&& y) {
  const auto prec = std::numeric_limits<std::remove_reference_t<T1>>::epsilon();
  return almost_eq(x, y, prec);
}

/**
 * @brief 在给定精度下判断 x 是否近似小于等于 y。
 * @param [in] x 第一个浮点数。
 * @param [in] y 第二个浮点数。
 * @param [in] prec 比较精度。
 * @return x <= y 或近似相等时返回 true。
 */
template <class T1, class T2, class T3>
bool almost_le(T1&& x, T2&& y, T3&& prec) {
  return x < y || almost_eq(x, y, prec);
}

/**
 * @brief 使用机器精度判断 x 是否近似小于等于 y。
 * @param [in] x 第一个浮点数。
 * @param [in] y 第二个浮点数。
 * @return x <= y 或近似相等时返回 true。
 */
template <class T1, class T2, class T3>
constexpr bool almost_le(T1&& x, T2&& y) {
  return x < y || almost_eq(x, y);
}

/**
 * @brief 在给定精度下判断 x 是否近似大于等于 y。
 * @param [in] x 第一个浮点数。
 * @param [in] y 第二个浮点数。
 * @param [in] prec 比较精度。
 * @return x >= y 或近似相等时返回 true。
 */
template <class T1, class T2, class T3>
constexpr bool almost_ge(T1&& x, T2&& y, T3&& prec) {
  return x > y || almost_eq(x, y, prec);
}

/**
 * @brief 使用机器精度判断 x 是否近似大于等于 y。
 * @param [in] x 第一个浮点数。
 * @param [in] y 第二个浮点数。
 * @return x >= y 或近似相等时返回 true。
 */
template <class T1, class T2>
constexpr bool almost_ge(T1&& x, T2&& y) {
  return x > y || almost_eq(x, y);
}
}  // namespace numerics