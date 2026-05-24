/**
 * @file PerformanceIndexEvaluator.hpp
 * @brief 根据动力学模式汇总 rollout metrics 为 PerformanceIndex。
 */
#pragma once

#include <array>

#include "PerformanceIndex.hpp"
#include "TrapezoidalIntegration.hpp"
#include "iLQRDescriptorTraits.hpp"

/**
 * @brief PerformanceIndex 汇总器。
 * @tparam Descriptor iLQR 描述类型。
 * @tparam DynamicsMode 动力学模式 tag。
 */
template <typename Descriptor, typename DynamicsMode>
struct PerformanceIndexEvaluator;

/** @brief 连续模式：按时间戳对 running metrics 做梯形积分。 */
template <typename Descriptor>
struct PerformanceIndexEvaluator<Descriptor, ContinuousDynamics> {
  using Traits = iLQRDescriptorTraits<Descriptor>;
  using Scalar = typename Traits::Scalar;

  static constexpr std::size_t PredictLength = Traits::PredictLength;

  using TimeTrajectory_t = typename Traits::TimeTrajectory_t;
  using PerformanceIndex_t = PerformanceIndex<Scalar>;
  using IntermediatePerformanceIndexTrajectory_t =
      std::array<PerformanceIndex_t, PredictLength>;

  static PerformanceIndex_t evaluate(
      const TimeTrajectory_t& timeTrajectory,
      const IntermediatePerformanceIndexTrajectory_t& intermediates,
      const PerformanceIndex_t& finalPerformanceIndex) {
    return trapezoidalIntegration(timeTrajectory, intermediates,
                                  finalPerformanceIndex);
  }
};

/** @brief 离散模式：直接累加 stage metrics 与终端 metrics。 */
template <typename Descriptor>
struct PerformanceIndexEvaluator<Descriptor, DiscreteDynamics> {
  using Traits = iLQRDescriptorTraits<Descriptor>;
  using Scalar = typename Traits::Scalar;

  static constexpr std::size_t PredictLength = Traits::PredictLength;

  using TimeTrajectory_t = typename Traits::TimeTrajectory_t;
  using PerformanceIndex_t = PerformanceIndex<Scalar>;
  using IntermediatePerformanceIndexTrajectory_t =
      std::array<PerformanceIndex_t, PredictLength>;

  static PerformanceIndex_t evaluate(
      const TimeTrajectory_t& timeTrajectory,
      const IntermediatePerformanceIndexTrajectory_t& intermediates,
      const PerformanceIndex_t& finalPerformanceIndex) {
    (void)timeTrajectory;
    PerformanceIndex_t result = finalPerformanceIndex;
    for (const auto& intermediate : intermediates) {
      result += intermediate;
    }
    return result;
  }
};

/** @brief 给定 descriptor 的 PerformanceIndex 汇总器类型。 */
template <typename Descriptor>
using PerformanceIndexEvaluator_t = PerformanceIndexEvaluator<
    Descriptor, typename iLQRDescriptorTraits<Descriptor>::DynamicsMode>;
