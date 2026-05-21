/******************************************************************************
Copyright (c) 2017, Farbod Farshidian. All rights reserved.

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
  }  // end of k loop

  return initialValue;
}