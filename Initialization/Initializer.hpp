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
 * @file Initializer.hpp
 * @brief 轨迹初始化器接口：在无控制器的时间段内，由 (time, state) 与 nextTime
 * 计算 input 与 nextState。
 */
#pragma once

#include "Types.hpp"

/**
 * @brief 求解器在无控制器可用时使用的初始化器接口；简单实现见
 * DefaultInitializer。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class Initializer {
 public:
  Initializer() = default;
  virtual ~Initializer() = default;

  /**
   * @brief 由当前时间与状态及下一时刻，计算当前输入与下一状态。
   * @param [in] time 当前时间。
   * @param [in] state 当前状态。
   * @param [in] nextTime 下一时刻（通常为 time + timeStep）。
   * @param [out] input 当前段的输入。
   * @param [out] nextState 下一时刻的状态。
   */
  virtual void compute(const Scalar time, const Vector<Scalar, XDim>& state,
                       const Scalar nextTime, Vector<Scalar, UDim>& input,
                       Vector<Scalar, XDim>& nextState) = 0;

 protected:
  /** @brief 拷贝构造（保护）。 */
  Initializer(const Initializer& rhs) = default;
};
