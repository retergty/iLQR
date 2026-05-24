/**
 * @file TrapezoidalIntegration.hpp
 * @brief 梯形积分：对时间-值轨迹做梯形法则积分，需 VALUE 支持加法与标量乘法。
 */
#pragma once
#include <stdint.h>

#include <array>

/**
 * @brief 对 (timeTrajectory, valueTrajectory) 做梯形积分，从 initialValue
 * 起累加。
 * @param [in] timeTrajectory 时间序列。
 * @param [in] valueTrajectory 值序列。
 * @param [in] initialValue 初始累加值。
 * @return 梯形积分结果。
 */
template <typename Scalar, typename VALUE, std::size_t TimeArrayLen,
          std::size_t ValueArrayLen>
VALUE trapezoidalIntegration(
    const std::array<Scalar, TimeArrayLen>& timeTrajectory,
    const std::array<VALUE, ValueArrayLen>& valueTrajectory,
    VALUE initialValue) {
  constexpr std::size_t ArrayLen = std::min(TimeArrayLen, ValueArrayLen);

  if constexpr (ArrayLen < 2) {
    return initialValue;
  }

  for (std::size_t k = 1; k < ArrayLen; k++) {
    VALUE temp = valueTrajectory[k - 1] + valueTrajectory[k];
    temp *= (0.5 * (timeTrajectory[k] - timeTrajectory[k - 1]));
    initialValue += temp;
  }  // k 循环结束。

  return initialValue;
}