#pragma once

#include <array>

#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "QpDiscreteTranscription.hpp"
#include "QpSolver.hpp"
#include "QpTrajectories.hpp"

namespace qp_solver {

/**
 * 围绕给定线性化轨迹求解带约束的离散时间线性二次控制问题。
 * 时间区间和离散步数由给定线性化轨迹的时间轨迹定义。
 *
 * @param optimalControlProblem: 最优控制问题定义。
 * @param nominalTrajectory: 用于构造线性二次近似的时间、状态和输入轨迹。
 * @param initialState: 预测区间起点状态。
 * @param intermediateMultipliers: 中间增广拉格朗日项的乘子。
 * 拉格朗日项。
 * @param finalMultipliers: 终端增广拉格朗日项的乘子。
 * @return 时间、状态和输入解。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                     Transcription::PredictLength>
solveLinearQuadraticOptimalControlProblem(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const TargetTrajectories<Scalar, Transcription>& targetTrajectory,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& nominalTrajectory,
    const Vector<Scalar, Transcription::XDim>& initialState,
    const std::array<
        MultiplierCollection<
            Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
        Transcription::PredictLength>& intermediateMultipliers,
    const MultiplierCollection<Scalar,
                               FinalStageConstraintLayout<ConstraintConfig>>&
        finalMultipliers) {
  constexpr int XDim = Transcription::XDim;
  constexpr int UDim = Transcription::UDim;
  constexpr size_t PredictLength = Transcription::PredictLength;

  // 近似
  const auto lqApproximation = getLinearQuadraticApproximation(
      optimalControlProblem, targetTrajectory, nominalTrajectory,
      intermediateMultipliers, finalMultipliers);

  // 求解更新步
  ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> deltaSolution;
  deltaSolution.timeTrajectory = nominalTrajectory.timeTrajectory;
  const Vector<Scalar, XDim> dx0 =
      initialState - nominalTrajectory.stateTrajectory.front();
  const auto deltaTrajectories =
      solveLinearQuadraticProblem(lqApproximation, dx0);
  const auto& stateDeltaTrajectory = deltaTrajectories.first;
  const auto& inputDeltaTrajectory = deltaTrajectories.second;

  for (size_t k = 0; k < PredictLength; ++k) {
    deltaSolution.stateTrajectory[k] = stateDeltaTrajectory[k];
    deltaSolution.inputTrajectory[k] = inputDeltaTrajectory[k];
  }
  deltaSolution.stateTrajectory[PredictLength] =
      stateDeltaTrajectory[PredictLength];

  // 执行完整步长：将更新量加入名义轨迹。
  return nominalTrajectory + deltaSolution;
}

/**
 * 用于未外部提供乘子轨迹的问题的重载。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                     Transcription::PredictLength>
solveLinearQuadraticOptimalControlProblem(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const TargetTrajectories<Scalar, Transcription>& targetTrajectory,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& nominalTrajectory,
    const Vector<Scalar, Transcription::XDim>& initialState) {
  std::array<MultiplierCollection<
                 Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
             Transcription::PredictLength>
      intermediateMultipliers{};
  MultiplierCollection<Scalar, FinalStageConstraintLayout<ConstraintConfig>>
      finalMultipliers{};
  return solveLinearQuadraticOptimalControlProblem(
      optimalControlProblem, targetTrajectory, nominalTrajectory, initialState,
      intermediateMultipliers, finalMultipliers);
}

template <typename Scalar, typename Transcription, typename ConstraintConfig>
ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                     Transcription::PredictLength>
solveLinearQuadraticOptimalControlProblem(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& nominalTrajectory,
    const Vector<Scalar, Transcription::XDim>& initialState) {
  std::array<MultiplierCollection<
                 Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
             Transcription::PredictLength>
      intermediateMultipliers{};
  MultiplierCollection<Scalar, FinalStageConstraintLayout<ConstraintConfig>>
      finalMultipliers{};
  const auto targetTrajectory =
      toTargetTrajectories<Scalar, Transcription>(nominalTrajectory);
  return solveLinearQuadraticOptimalControlProblem(
      optimalControlProblem, targetTrajectory, nominalTrajectory, initialState,
      intermediateMultipliers, finalMultipliers);
}

}  // namespace qp_solver
