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
 * @file QuadraticPenalty.hpp
 * @brief 二次惩罚：等式约束 h=0 的增广拉格朗日惩罚 p(h,λ)=λh+(μ/2)h²，乘子更新 λ_{k+1}=λ_k-α*h。
 */
#pragma once

#include "AugmentedPenaltyBase.hpp"

 /**
  * @brief 单等式约束的二次惩罚实现：惩罚项 p = λh + (μ/2)h²，μ 为尺度参数。
  * @tparam Scalar 标量类型。
  */
template<typename Scalar>
class QuadraticPenalty final : public AugmentedPenaltyBase<Scalar>
{
public:
  /**
   * Configuration object for the quadratic penalty.
   * scale: scaling factor, see class description
   * stepSize: step-length parameter, see class description
   */
  struct Config {
    Config() : Config(100.0, 0.0) {}
    Config(const Scalar scaleParam, const Scalar stepSizeParam) : scale(scaleParam), stepSize(stepSizeParam) {}
    Scalar scale;
    Scalar stepSize;
  };

  /** Constructor */
  explicit QuadraticPenalty(const Config& config) : config_(config) {}

  ~QuadraticPenalty() override = default;

  Scalar getValue(const Scalar t, const Scalar l, const Scalar h) const override {
    (void)t;
    return -l * h + 0.5 * config_.scale * h * h;
  }
  Scalar getDerivative(const Scalar t, const Scalar l, const Scalar h) const override
  {
    (void)t;
    return -l + config_.scale * h;
  }
  Scalar getSecondDerivative(const Scalar t, const Scalar l, const Scalar h) const override {
    (void)t;
    (void)l;
    (void)h;
    return config_.scale;
  }

  Scalar updateMultiplier(const Scalar t, const Scalar l, const Scalar h) const override {
    (void)t;
    return l - config_.stepSize * config_.scale * h;
  }
  Scalar initializeMultiplier() const override { return 0.0; }

private:
  QuadraticPenalty(const QuadraticPenalty& other) = default;

  const Config config_;
};
