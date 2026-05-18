/******************************************************************************
Copyright (c) 2020, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

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
constexpr bool almost_eq(T1 &&x, T2 &&y, T3 &&prec) {
  static_assert(
      std::is_floating_point<typename std::remove_reference<T1>::type>::value,
      "First argument is not floating point!");
  static_assert(
      std::is_floating_point<typename std::remove_reference<T2>::type>::value,
      "Second argument is not floating point!");
  static_assert(
      std::is_floating_point<typename std::remove_reference<T3>::type>::value,
      "prec is not floating point!");
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision unless the result is subnormal
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
template <class T1, class T2> constexpr bool almost_eq(T1 &&x, T2 &&y) {
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
bool almost_le(T1 &&x, T2 &&y, T3 &&prec) {
  return x < y || almost_eq(x, y, prec);
}

/**
 * @brief 使用机器精度判断 x 是否近似小于等于 y。
 * @param [in] x 第一个浮点数。
 * @param [in] y 第二个浮点数。
 * @return x <= y 或近似相等时返回 true。
 */
template <class T1, class T2, class T3>
constexpr bool almost_le(T1 &&x, T2 &&y) {
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
constexpr bool almost_ge(T1 &&x, T2 &&y, T3 &&prec) {
  return x > y || almost_eq(x, y, prec);
}

/**
 * @brief 使用机器精度判断 x 是否近似大于等于 y。
 * @param [in] x 第一个浮点数。
 * @param [in] y 第二个浮点数。
 * @return x >= y 或近似相等时返回 true。
 */
template <class T1, class T2> constexpr bool almost_ge(T1 &&x, T2 &&y) {
  return x > y || almost_eq(x, y);
}
} // namespace numerics