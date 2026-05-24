/**
 * @file LinearInterpolation.hpp
 * @brief 线性插值：在固定步长时间数组上计算区间、插值系数及轨迹插值。
 */
#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <type_traits>

namespace LinearInterpolation {
/** @brief 在严格等间隔时间数组中用 O(1) 算术定位查询时间所在区间。 */
template <typename Scalar, size_t ArrayLength>
std::pair<int, Scalar> timeSegment(
    Scalar enquiryTime, const std::array<Scalar, ArrayLength>& timeArray) {
  // 边界情况（未设置时间或只有一个时间元素）。
  if constexpr (ArrayLength <= 1) {
    return {0, Scalar(1.0)};
  } else {
    const Scalar dt = timeArray[1] - timeArray[0];
    assert(dt > Scalar(0));

    constexpr int lastInterval = static_cast<int>(ArrayLength) - 1;
    constexpr int maxIndex = lastInterval - 1;

    const Scalar position = (enquiryTime - timeArray.front()) / dt;
    int index = static_cast<int>(std::floor(position));
    index = std::clamp(index, 0, maxIndex);

    Scalar alpha = Scalar(1) - (position - static_cast<Scalar>(index));
    alpha = std::clamp(alpha, Scalar(0), Scalar(1));
    return {index, alpha};
  }
}

template <typename Data, size_t ArrayLength>
const Data& stdAccessFun(const std::array<Data, ArrayLength>& arr, size_t ind) {
  return arr[ind];
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar, typename Data, size_t ArrayLength, class AccessFun>
auto interpolate(const std::pair<int, Scalar>& indexAlpha,
                 const std::array<Data, ArrayLength>& dataArray,
                 AccessFun accessFun)
    -> std::decay_t<std::invoke_result_t<
        AccessFun, const std::array<Data, ArrayLength>&, size_t>> {
  static_assert(ArrayLength > 0);

  if constexpr (ArrayLength > 1) {
    // 普通插值情况。
    int index = indexAlpha.first;
    Scalar alpha = indexAlpha.second;
    const auto& lhs = accessFun(dataArray, index);
    const auto& rhs = accessFun(dataArray, index + 1);

    return alpha * lhs + (Scalar(1.0) - alpha) * rhs;
  } else {  // dataArray.size() == 1
            // 时间向量只有 1 个元素 -> 常值函数。
    return accessFun(dataArray, 0);
  }
}

template <typename Scalar, typename Data, size_t ArrayLength, class AccessFun>
auto interpolate(const Scalar enquiryTime,
                 const std::array<Scalar, ArrayLength>& timeArray,
                 const std::array<Data, ArrayLength>& dataArray,
                 AccessFun accessFun)
    -> std::decay_t<std::invoke_result_t<
        AccessFun, const std::array<Data, ArrayLength>&, size_t>> {
  return interpolate(timeSegment(enquiryTime, timeArray), dataArray, accessFun);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar, typename Data, size_t ArrayLength>
Data interpolate(const std::pair<int, Scalar>& indexAlpha,
                 const std::array<Data, ArrayLength>& dataArray) {
  return interpolate(indexAlpha, dataArray, stdAccessFun<Data, ArrayLength>);
}

template <typename Scalar, typename TimeType, typename Data, size_t ArrayLength>
Data interpolate(const Scalar enquiryTime,
                 const std::array<TimeType, ArrayLength>& timeArray,
                 const std::array<Data, ArrayLength>& dataArray) {
  return interpolate(enquiryTime, timeArray, dataArray,
                     stdAccessFun<Data, ArrayLength>);
}

}  // namespace LinearInterpolation