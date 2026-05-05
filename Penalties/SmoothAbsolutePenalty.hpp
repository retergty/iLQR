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
 * @file SmoothAbsolutePenalty.hpp
 * @brief 光滑绝对值惩罚：用于等式约束 h=0，p(h)=μ√(h²+δ²)，乘子更新 λ_{k+1}=λ_k-α*h。
 */
#pragma once

#include "AugmentedPenaltyBase.hpp"

/**
 * @brief 单等式约束的光滑绝对值惩罚实现：p(h)=μ√(h²+δ²)，μ 为尺度、δ 为松弛参数。
 * @tparam Scalar 标量类型。
 */
template<typename Scalar>
class SmoothAbsolutePenalty final : public AugmentedPenaltyBase {
public:
  /**
   * Configuration object for the smooth absolute penalty.
   * scale: scaling factor, see class description
   * relaxation: relaxation parameter, see class description
   * stepLenght: step-length parameter, see class description
   */
  struct Config {
    Config() : Config(100.0, 1e-2, 0.0) {}
    Config(const Scalar scaleParam, const Scalar relaxationParam, const Scalar stepSizeParam)
      : scale(scaleParam), relaxation(relaxationParam), stepSize(stepSizeParam) {
    }
    Scalar scale;
    Scalar relaxation;
    Scalar stepSize;
  };

  /** Constructor */
  explicit SmoothAbsolutePenalty(const Config &config) : config_(config) {}

  ~SmoothAbsolutePenalty() override = default;

  Scalar getValue(const Scalar t, const Scalar l, const Scalar h) const override {
    return -l * h + config_.scale * sqrt(h * h + config_.relaxation * config_.relaxation);
  }
  Scalar getDerivative(const Scalar t, const Scalar l, const Scalar h) const override {
    return -l + config_.scale * h / sqrt(h * h + config_.relaxation * config_.relaxation);
  }
  Scalar getSecondDerivative(const Scalar t, const Scalar l, const Scalar h) const override {
    const const Scalar deltaSquare = config_.relaxation * config_.relaxation;
    return config_.scale * deltaSquare / pow(h * h + deltaSquare, 1.5);
  }

  Scalar updateMultiplier(const Scalar t, const Scalar l, const Scalar h) const override { return l - config_.stepSize * config_.scale * h; }
  Scalar initializeMultiplier() const override { return 0.0; }

private:
  SmoothAbsolutePenalty(const SmoothAbsolutePenalty& other) = default;

  const Config config_;
};