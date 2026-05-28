/**
 * @file ProblemMetrics.hpp
 * @brief 整条轨迹的问题指标：各中间时刻与终端时刻的 Metrics 汇总。
 */
#pragma once

#include "ModelData/Metrics.hpp"
#include "iLQR/iLQRDescriptor.hpp"

/**
 * @brief 整条 rollout 的问题指标容器：终端一点 Metrics + 中间 PredictLength
 * 个点的 Metrics。
 * @tparam Scalar 标量类型。
 * @tparam Transcription 轨迹配置，提供 XDim/UDim/PredictLength。
 * @tparam ConstraintConfig 约束配置。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
struct ProblemMetrics {
  static constexpr std::size_t PredictLength = Transcription::PredictLength;

  using IntermediateMetrics_t =
      Metrics<Scalar, typename Transcription::Dims,
              IntermediateStageConstraintLayout<ConstraintConfig>>;
  using FinalMetrics_t = Metrics<Scalar, typename Transcription::Dims,
                                 FinalStageConstraintLayout<ConstraintConfig>>;

  /** @brief 终端时刻的 Metrics。 */
  FinalMetrics_t final;
  /** @brief 各中间时刻的 Metrics，长度为 PredictLength。 */
  std::array<IntermediateMetrics_t, PredictLength> intermediates;

  /** @brief 与另一 ProblemMetrics 交换内容。 */
  void swap(ProblemMetrics& other) {
    final.swap(other.final);
    intermediates.swap(other.intermediates);
  }

  /** @brief 清空/重置各 Metrics（用于重新初始化）。 */
  void clear() {
    final.clear();
    for (size_t i = 0; i < intermediates.size(); ++i) intermediates[i].clear();
  }
};