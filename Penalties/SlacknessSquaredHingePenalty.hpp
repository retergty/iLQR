/**
 * @file SlacknessSquaredHingePenalty.hpp
 * @brief 松弛平方铰链惩罚（PHR）：不等式 h≥0 的增广拉格朗日惩罚，引入松弛 s 后
 * p=(1/2ρ)(max(0,λ-ρh)²-λ²)。
 */
#pragma once

#include "Penalties/AugmentedPenaltyBase.hpp"

/**
 * @brief 单不等式约束的 PHR 惩罚实现：乘子更新为 λ_{k+1}=max(λ_k-α*h,
 * (1-α/ρ)λ_k)。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
class SlacknessSquaredHingePenalty final : public AugmentedPenaltyBase<Scalar> {
 public:
  /**
   * 平方铰链惩罚的配置对象。
   * @param [in] scale: 缩放因子，在类说明中记为
   * \pho。
   * @param [in] stepSize: 更新拉格朗日乘子的步长，在
   * 类说明中记为 \alpha。
   */
  struct Config {
    Config() : Config(Scalar(10.0), Scalar(1.0)) {}
    Config(const Scalar scaleParam, const Scalar stepSizeParam)
        : scale(scaleParam), stepSize(stepSizeParam) {}
    Scalar scale;
    Scalar stepSize;
  };

  /** 构造函数 */
  SlacknessSquaredHingePenalty(const Config& config) : config_(config) {}

  ~SlacknessSquaredHingePenalty() override = default;

  Scalar getValue(const Scalar t, const Scalar l,
                  const Scalar h) const override {
    (void)t;
    return (h < l / config_.scale)
               ? (-l * h + Scalar(0.5) * config_.scale * h * h)
               : (-Scalar(0.5) * l * l / config_.scale);
  }
  Scalar getDerivative(const Scalar t, const Scalar l,
                       const Scalar h) const override {
    (void)t;
    return (h < l / config_.scale) ? (-l + config_.scale * h) : Scalar(0.0);
  }
  Scalar getSecondDerivative(const Scalar t, const Scalar l,
                             const Scalar h) const override {
    (void)t;
    return (h < l / config_.scale) ? config_.scale : Scalar(0.0);
  }

  Scalar updateMultiplier(const Scalar t, const Scalar l,
                          const Scalar h) const override {
    (void)t;
    return std::max(Scalar(0.0),
                    std::max(l - config_.stepSize * config_.scale * h,
                             (Scalar(1.0) - config_.stepSize) * l));
  }
  Scalar initializeMultiplier() const override { return Scalar(0.0); }

 private:
  SlacknessSquaredHingePenalty(const SlacknessSquaredHingePenalty& other) =
      default;

  const Config config_;
};