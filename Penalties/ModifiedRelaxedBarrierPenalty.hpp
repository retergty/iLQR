/**
 * @file ModifiedRelaxedBarrierPenalty.hpp
 * @brief 修正松弛对数障碍惩罚：不等式 h≥0，p=(λ²/ρ)ψ(ρh/λ)，ψ 为平移二次松弛
 * log barrier，定义域 x>-1。
 */
#pragma once

#include <math.h>

#include "Penalties/AugmentedPenaltyBase.hpp"

/**
 * @brief 单不等式约束的修正松弛障碍惩罚实现：乘子更新
 * λ_{k+1}=-α*λ_k*ψ'(ρh/λ_k)。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
class ModifiedRelaxedBarrierPenalty final
    : public AugmentedPenaltyBase<Scalar> {
 public:
  /**
   * modified relaxed barrier 惩罚的配置对象。
   * scale：缩放因子，参见类说明。
   * relaxation：松弛参数，参见类说明。
   * stepLenght：步长参数，参见类说明。
   */
  struct Config {
    Config() : Config(Scalar(10.0), Scalar(0.0), Scalar(1.0)) {}
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
  explicit ModifiedRelaxedBarrierPenalty(const Config& config)
      : config_(config), quadCoeff_(config) {}

  ~ModifiedRelaxedBarrierPenalty() override = default;

  Scalar getValue(const Scalar t, const Scalar l,
                  const Scalar h) const override {
    (void)t;
    const Scalar v = vFunc(l, h);
    if (v > config_.relaxation) {
      return -wFunc(l) * log(Scalar(1.0) + v);
    } else {
      const Scalar vDelta = v - config_.relaxation;
      return wFunc(l) *
             (static_cast<Scalar>(0.5) * quadCoeff_.c2 * vDelta * vDelta +
              quadCoeff_.c1 * vDelta + quadCoeff_.c0);
    }
  }

  Scalar getDerivative(const Scalar t, const Scalar l,
                       const Scalar h) const override {
    (void)t;
    const Scalar v = vFunc(l, h);
    if (v > config_.relaxation) {
      return -wFunc(l) / (Scalar(1.0) + v) * dvdhFunc(l);
    } else {
      return wFunc(l) *
             (quadCoeff_.c2 * (v - config_.relaxation) + quadCoeff_.c1) *
             dvdhFunc(l);
    }
  }

  Scalar getSecondDerivative(const Scalar t, const Scalar l,
                             const Scalar h) const override {
    (void)t;
    const Scalar v = vFunc(l, h);
    const Scalar dvdh = dvdhFunc(l);
    if (v > config_.relaxation) {
      return wFunc(l) / ((Scalar(1.0) + v) * (Scalar(1.0) + v)) * dvdh * dvdh;
    } else {
      return wFunc(l) * quadCoeff_.c2 * dvdh * dvdh;
    }
  }

  Scalar updateMultiplier(const Scalar t, const Scalar l,
                          const Scalar h) const override {
    (void)t;
    const Scalar v = vFunc(l, h);
    constexpr Scalar lambdaMin = Scalar(1e-4);
    if (v > config_.relaxation) {
      return std::max(lambdaMin, wFunc(l) * dvdhFunc(l) / (Scalar(1.0) + v));
    } else {
      return std::max(
          lambdaMin,
          config_.stepSize * wFunc(l) *
              (-quadCoeff_.c2 * (v - config_.relaxation) - quadCoeff_.c1) *
              dvdhFunc(l));
    }
  }

  Scalar initializeMultiplier() const override { return Scalar(1.0); }

 private:
  ModifiedRelaxedBarrierPenalty(const ModifiedRelaxedBarrierPenalty& other) =
      default;

  Scalar wFunc(const Scalar l) const { return l * l / config_.scale; }
  Scalar dvdhFunc(const Scalar l) const { return config_.scale / l; }
  Scalar vFunc(const Scalar l, const Scalar h) const {
    return config_.scale * h / l;
  }

  struct QuadCoeff {
    QuadCoeff(const Config& config) {
      c2 = Scalar(1.0) / std::pow(Scalar(1.0) + config.relaxation, Scalar(2.0));
      c1 = -Scalar(1.0) / (Scalar(1.0) + config.relaxation);
      c0 = -log(Scalar(1.0) + config.relaxation);
    }

    Scalar c2 = Scalar(0.0), c1 = Scalar(0.0), c0 = Scalar(0.0);
  };

  const Config config_;
  const QuadCoeff quadCoeff_;
};
