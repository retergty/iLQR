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
 * @file OperatingPoints.hpp
 * @brief 工作点初始化器：用固定状态/输入工作点生成轨迹，不依赖动力学积分。
 */
#pragma once

#include "Initializer.hpp"

/**
 * @brief 基于工作点的初始化器：输出恒为给定状态工作点与输入工作点，nextState
 * 为状态工作点。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class OperatingPoints final : public Initializer<Scalar, XDim, UDim> {
 public:
  /**
   * @brief 用状态工作点与输入工作点构造。
   * @param [in] stateOperatingPoint 状态工作点。
   * @param [in] inputOperatingPoint 输入工作点。
   */
  OperatingPoints(const Vector<Scalar, XDim>& stateOperatingPoint,
                  const Vector<Scalar, UDim>& inputOperatingPoint)
      : stateTrajectory_(stateOperatingPoint),
        inputTrajectory_(inputOperatingPoint) {}

  /** @brief 析构函数。 */
  ~OperatingPoints() override = default;

  /** @brief 将 input 设为输入工作点，nextState 设为状态工作点。 */
  void compute(const Scalar time, const Vector<Scalar, XDim>& state,
               const Scalar nextTime, Vector<Scalar, UDim>& input,
               Vector<Scalar, XDim>& nextState) override {
    (void)time;
    (void)state;
    (void)nextTime;
    input = inputTrajectory_;
    nextState = stateTrajectory_;
  }

 private:
  /** @brief 拷贝构造（保护）。 */
  OperatingPoints(const OperatingPoints& other) = default;

  /** @brief 状态工作点。 */
  const Vector<Scalar, XDim> stateTrajectory_;
  /** @brief 输入工作点。 */
  const Vector<Scalar, UDim> inputTrajectory_;
};
