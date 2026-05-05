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
 * @file DefaultInitializer.hpp
 * @brief 默认初始化器：输入置零，下一状态等于当前状态（保持不动）。
 */
#pragma once

#include "Initializer.hpp"

/**
 * @brief 默认初始化器实现：input 置零，nextState = state。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template<typename Scalar, int XDim, int UDim>
class DefaultInitializer final : public Initializer<Scalar, XDim, UDim>
{
public:
  /** @brief 默认构造。 */
  explicit DefaultInitializer() = default;

  ~DefaultInitializer() override = default;

  /** @brief 将 input 置零，nextState 设为当前 state。 */
  void compute(const Scalar time, const Vector<Scalar, XDim>& state, const Scalar nextTime, Vector<Scalar, UDim>& input, Vector<Scalar, XDim>& nextState) override {
    (void)time;
    (void)nextTime;
    (void)state;
    
    input.setZero();
    nextState = state;
  }

protected:
  DefaultInitializer(const DefaultInitializer& rhs) = default;
};
