/**
 * @file DDPData.hpp
 * @brief DDP 求解器用的原始/对偶数据容器：PrimalDataContainer 与
 * DualDataContainer。
 */
#pragma once

#include "DualSolution.hpp"
#include "ModelData.hpp"
#include "PrimalSolution.hpp"
#include "ProblemMetrics.hpp"
#include "RiccatiModification.hpp"

/**
 * @brief 原始数据容器：存放一次 rollout 的轨迹、控制器、中间/终端模型数据与
 * 问题指标。
 * @note 各轨迹与 modelData 应来自同一控制器的
 * rollout；用外部控制器初始化时需随后调用 run 以补齐数据。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
struct PrimalDataContainer {
  static constexpr int XDim = Transcription::XDim;
  static constexpr int UDim = Transcription::UDim;
  static constexpr std::size_t PredictLength = Transcription::PredictLength;

  using PrimalSolution_t = PrimalSolution<Scalar, Transcription>;
  using ProblemMetrics_t =
      ProblemMetrics<Scalar, Transcription, ConstraintConfig>;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;

  // 原始解
  PrimalSolution_t primalSolution;
  // rollout 的代价、软约束和约束值。
  ProblemMetrics_t problemMetrics;
  // 终端时刻模型数据。
  ModelData_t modelDataFinalTime;
  // 中间节点模型数据轨迹。
  std::array<ModelData_t, PredictLength> modelDataTrajectory;

  /** @brief 与另一容器交换内容。 */
  void swap(PrimalDataContainer& other) {
    primalSolution.swap(other.primalSolution);
    problemMetrics.swap(other.problemMetrics);
    std::swap(modelDataFinalTime, other.modelDataFinalTime);
    modelDataTrajectory.swap(other.modelDataTrajectory);
  }

  /** @brief 清空 primal 解与 problem metrics；modelDataTrajectory
   * 在下次近似时被覆盖。 */
  void clear() {
    primalSolution.clear();
    problemMetrics.clear();
    // std::array 没有 clear()；modelDataTrajectory 会在下一次
    // approximateOptimalControlProblem() 中被覆盖。
  }
};

/**
 * @brief 对偶数据容器：存放对偶解、Riccati 修正轨迹与 value
 * function 轨迹。
 * @note valueFunctionTrajectory 由模型数据与 riccatiModification 经
 * Riccati 递推得到。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
struct DualDataContainer {
  static constexpr int XDim = Transcription::XDim;
  static constexpr int UDim = Transcription::UDim;
  static constexpr std::size_t PredictLength = Transcription::PredictLength;

  using DualSolution_t =
      DualSolution<Scalar, typename Transcription::Horizon, ConstraintConfig>;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;
  using RiccatiModification_t = RiccatiModification<Scalar, XDim, UDim>;
  using ValueFunctionQuadraticApproximation_t =
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>;
  // 对偶解
  DualSolution_t dualSolution;
  // Riccati 修正。
  std::array<RiccatiModification_t, PredictLength + 1>
      riccatiModificationTrajectory;
  // Riccati 解系数。
  std::array<ValueFunctionQuadraticApproximation_t, PredictLength + 1>
      valueFunctionTrajectory;

  /** @brief 与另一容器交换内容。 */
  void swap(DualDataContainer& other) {
    dualSolution.swap(other.dualSolution);
    riccatiModificationTrajectory.swap(other.riccatiModificationTrajectory);
    valueFunctionTrajectory.swap(other.valueFunctionTrajectory);
  }
};
