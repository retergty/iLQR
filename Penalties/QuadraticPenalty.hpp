#pragma once

#include "AugmentedPenaltyBase.hpp"

/**
 * @brief 等式约束 \f$ h = 0 \f$ 的二次增广惩罚项。
 *
 * 该类实现单个标量等式约束的增广拉格朗日惩罚函数：
 *
 * \f[
 *   p(h, \lambda) = -\lambda h + \frac{\mu}{2} h^2
 * \f]
 *
 * 其中 \f$ h \f$ 是约束值，\f$ \lambda \f$ 是拉格朗日乘子，
 * \f$ \mu \f$ 对应配置中的 `scale`，用于控制二次惩罚强度。
 *
 * 外层 `StateAugmentedLagrangian` / `StateInputAugmentedLagrangian` 会再用
 * `Multiplier::penalty` 对该惩罚值整体缩放。当前实现中的实际代价贡献为：
 *
 * \f[
 *   L_A = \rho \, p(h, \lambda)
 * \f]
 *
 * 其中 \f$ \rho \f$ 是 `Multiplier::penalty`。
 *
 * 乘子更新规则为：
 *
 * \f[
 *   \lambda_{k+1} = \lambda_k - \eta \mu h_{k+1}
 * \f]
 *
 * 其中 \f$ \eta \f$ 对应配置中的 `stepSize`。因此 `stepSize = 1`
 * 时，更新步长等于二次惩罚系数 \f$ \mu \f$。
 */
template <typename Scalar>
class QuadraticPenalty final : public AugmentedPenaltyBase<Scalar> {
 public:
  /**
   * Configuration object for the quadratic penalty.
   * scale: scaling factor, see class description
   * stepSize: step-length parameter, see class description
   */
  struct Config {
    Config() : Config(100.0, 0.0) {}
    Config(const Scalar scaleParam, const Scalar stepSizeParam)
        : scale(scaleParam), stepSize(stepSizeParam) {}
    Scalar scale;
    Scalar stepSize;
  };

  /** Constructor */
  explicit QuadraticPenalty(const Config& config) : config_(config) {}

  ~QuadraticPenalty() override = default;

  Scalar getValue(const Scalar t, const Scalar l,
                  const Scalar h) const override {
    (void)t;
    return -l * h + 0.5 * config_.scale * h * h;
  }
  Scalar getDerivative(const Scalar t, const Scalar l,
                       const Scalar h) const override {
    (void)t;
    return -l + config_.scale * h;
  }
  Scalar getSecondDerivative(const Scalar t, const Scalar l,
                             const Scalar h) const override {
    (void)t;
    (void)l;
    (void)h;
    return config_.scale;
  }

  Scalar updateMultiplier(const Scalar t, const Scalar l,
                          const Scalar h) const override {
    (void)t;
    return l - config_.stepSize * config_.scale * h;
  }
  Scalar initializeMultiplier() const override { return 0.0; }

 private:
  QuadraticPenalty(const QuadraticPenalty& other) = default;

  const Config config_;
};
