/**
 * @file DDPData.hpp
 * @brief DDP 求解器用的 LQ 近似与 Riccati 递推数据。
 */
#pragma once

#include <array>
#include <utility>

#include "Approximation/QuadraticApproximation.hpp"
#include "ModelData/ModelData.hpp"
#include "RiccatiEquations/RiccatiModification.hpp"

/**
 * @brief LQ 近似数据：存放名义轨迹对应的中间/终端模型数据。
 * @note 原始解和 rollout 指标由求解器单独持有，本容器只保存围绕名义轨迹计算出的
 * LQ 模型近似。
 */
template <typename Scalar, typename Transcription>
struct LQApproximationData {
  static constexpr int XDim = Transcription::XDim;
  static constexpr int UDim = Transcription::UDim;
  static constexpr std::size_t PredictLength = Transcription::PredictLength;

  using ModelData_t = ModelData<Scalar, XDim, UDim>;

  // 终端时刻模型数据。
  ModelData_t modelDataFinalTime;
  // 中间节点模型数据轨迹。
  std::array<ModelData_t, PredictLength> modelDataTrajectory;

  /** @brief 与另一容器交换内容。 */
  void swap(LQApproximationData& other) {
    std::swap(modelDataFinalTime, other.modelDataFinalTime);
    modelDataTrajectory.swap(other.modelDataTrajectory);
  }
};

/**
 * @brief Riccati 递推数据：存放 Riccati 修正轨迹与 value function 轨迹。
 * @note valueFunctionTrajectory 由模型数据与 riccatiModification 经
 * Riccati 递推得到。
 */
template <typename Scalar, typename Transcription>
struct RiccatiData {
  static constexpr int XDim = Transcription::XDim;
  static constexpr int UDim = Transcription::UDim;
  static constexpr std::size_t PredictLength = Transcription::PredictLength;

  using RiccatiModification_t = RiccatiModification<Scalar, XDim, UDim>;
  using ValueFunctionQuadraticApproximation_t =
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>;
  // Riccati 修正。
  std::array<RiccatiModification_t, PredictLength>
      riccatiModificationTrajectory;
  // Riccati 解系数。
  std::array<ValueFunctionQuadraticApproximation_t, PredictLength + 1>
      valueFunctionTrajectory;

  /** @brief 与另一容器交换内容。 */
  void swap(RiccatiData& other) {
    riccatiModificationTrajectory.swap(other.riccatiModificationTrajectory);
    valueFunctionTrajectory.swap(other.valueFunctionTrajectory);
  }
};
