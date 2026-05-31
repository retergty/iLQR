/**
 * @file SmoothAbsolutePenalty.hpp
 * @brief 光滑绝对值惩罚：用于等式约束 h=0，p(h)=μ√(h²+δ²)，乘子更新
 * λ_{k+1}=λ_k-α*h。
 */
#pragma once

#include "Penalties/AugmentedPenaltyBase.hpp"

/**
 * @brief 单等式约束的光滑绝对值惩罚实现：p(h)=μ√(h²+δ²)，μ 为尺度、δ
 * 为松弛参数。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
class SmoothAbsolutePenalty final : public AugmentedPenaltyBase<Scalar> {
 public:
  /**
   * 光滑绝对值惩罚的配置对象。
   * scale：缩放因子，参见类说明。
   * relaxation：松弛参数，参见类说明。
   * stepLenght：步长参数，参见类说明。
   */
  struct Config {
    Config() : Config(Scalar(100.0), Scalar(1e-2), Scalar(0.0)) {}
    Config(const Scalar scaleParam, const Scalar relaxationParam,
           const Scalar stepSizeParam)
        : scale(scaleParam),
          relaxation(relaxationParam),
          stepSize(stepSizeParam) {}
    Scalar scale;
    Scalar relaxation;
    Scalar stepSize;
  };

  /** 构造函数 */
  explicit SmoothAbsolutePenalty(const Config& config) : config_(config) {}

  ~SmoothAbsolutePenalty() override = default;

  Scalar getValue(const Scalar t, const Scalar l,
                  const Scalar h) const override {
    (void)t;
    return -l * h + config_.scale *
                        sqrt(h * h + config_.relaxation * config_.relaxation);
  }
  Scalar getDerivative(const Scalar t, const Scalar l,
                       const Scalar h) const override {
    (void)t;
    return -l + config_.scale * h /
                    sqrt(h * h + config_.relaxation * config_.relaxation);
  }
  Scalar getSecondDerivative(const Scalar t, const Scalar l,
                             const Scalar h) const override {
    (void)t;
    (void)l;
    const Scalar deltaSquare = config_.relaxation * config_.relaxation;
    return config_.scale * deltaSquare / pow(h * h + deltaSquare, Scalar(1.5));
  }

  Scalar updateMultiplier(const Scalar t, const Scalar l,
                          const Scalar h) const override {
    (void)t;
    return l - config_.stepSize * config_.scale * h;
  }
  Scalar initializeMultiplier() const override { return Scalar(0.0); }

 private:
  SmoothAbsolutePenalty(const SmoothAbsolutePenalty& other) = default;

  const Config config_;
};