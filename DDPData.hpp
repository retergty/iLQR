/******************************************************************************
Copyright (c) 2021, Farbod Farshidian. All rights reserved.

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
 * problem metrics。
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

  // Primal solution
  PrimalSolution_t primalSolution;
  // cost, soft constraints and constraints values of the rollout
  ProblemMetrics_t problemMetrics;
  // final time model data
  ModelData_t modelDataFinalTime;
  // intermediate model data trajectory
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
    // std::array has no clear(); modelDataTrajectory is overwritten on next
    // approximateOptimalControlProblem()
  }
};

/**
 * @brief 对偶数据容器：存放对偶解、投影模型轨迹、Riccati 修正轨迹与 value
 * function 轨迹。
 * @note valueFunctionTrajectory 由 (projectedModelData, riccatiModification) 经
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
  // Dual solution
  DualSolution_t dualSolution;
  // projected model data trajectory
  std::array<ModelData_t, PredictLength + 1> projectedModelDataTrajectory;
  // Riccati modification
  std::array<RiccatiModification_t, PredictLength + 1>
      riccatiModificationTrajectory;
  // Riccati solution coefficients
  std::array<ValueFunctionQuadraticApproximation_t, PredictLength + 1>
      valueFunctionTrajectory;

  /** @brief 与另一容器交换内容。 */
  void swap(DualDataContainer& other) {
    dualSolution.swap(other.dualSolution);
    projectedModelDataTrajectory.swap(other.projectedModelDataTrajectory);
    riccatiModificationTrajectory.swap(other.riccatiModificationTrajectory);
    valueFunctionTrajectory.swap(other.valueFunctionTrajectory);
  }
};
