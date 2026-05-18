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
 * @file SlacknessSquaredHingePenalty.hpp
 * @brief 松弛平方铰链惩罚（PHR）：不等式 h≥0 的增广拉格朗日惩罚，引入松弛 s 后
 * p=(1/2ρ)(max(0,λ-ρh)²-λ²)。
 */
#pragma once

#include "AugmentedPenaltyBase.hpp"

/**
 * @brief 单不等式约束的 PHR 惩罚实现：乘子更新为 λ_{k+1}=max(λ_k-α*h,
 * (1-α/ρ)λ_k)。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
class SlacknessSquaredHingePenalty final : public AugmentedPenaltyBase<Scalar> {
 public:
  /**
   * Configuration object for the squared hinge penalty.
   * @param [in] scale : scaling factor. In the class description, it is
   * referred to as \pho.
   * @param [in] stepSize: step-size for updating Lagrange multiplier. In the
   * class description, it is referred to as \alpha.
   */
  struct Config {
    Config() : Config(10.0, 1.0) {}
    Config(const Scalar scaleParam, const Scalar stepSizeParam)
        : scale(scaleParam), stepSize(stepSizeParam) {}
    Scalar scale;
    Scalar stepSize;
  };

  /** Constructor */
  SlacknessSquaredHingePenalty(const Config& config) : config_(config) {}

  ~SlacknessSquaredHingePenalty() override = default;

  Scalar getValue(const Scalar t, const Scalar l,
                  const Scalar h) const override {
    (void)t;
    return (h < l / config_.scale) ? (-l * h + 0.5 * config_.scale * h * h)
                                   : (-0.5 * l * l / config_.scale);
  }
  Scalar getDerivative(const Scalar t, const Scalar l,
                       const Scalar h) const override {
    (void)t;
    return (h < l / config_.scale) ? (-l + config_.scale * h) : 0.0;
  }
  Scalar getSecondDerivative(const Scalar t, const Scalar l,
                             const Scalar h) const override {
    (void)t;
    return (h < l / config_.scale) ? config_.scale : 0.0;
  }

  Scalar updateMultiplier(const Scalar t, const Scalar l,
                          const Scalar h) const override {
    (void)t;
    return std::max(0.0, std::max(l - config_.stepSize * config_.scale * h,
                                  (1.0 - config_.stepSize) * l));
  }
  Scalar initializeMultiplier() const override { return 0.0; }

 private:
  SlacknessSquaredHingePenalty(const SlacknessSquaredHingePenalty& other) =
      default;

  const Config config_;
};