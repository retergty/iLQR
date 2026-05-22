/**
 * @file OptimalControlProblemHelperFunction.hpp
 * @brief 最优控制问题辅助函数：初始化乘子、对偶解及近似等。
 */
#pragma once
#include "DualSolution.hpp"
#include "Metrics.hpp"
#include "OptimalControlProblem.hpp"
#include "PrimalSolution.hpp"
#include "ProblemMetrics.hpp"

/**
 * @brief 初始化终端时刻的乘子集合（等式与不等式拉格朗日）。
 * @param [in] ocp 最优控制问题。
 * @param [in] time 终端时间。
 * @param [out] multiplierCollection 待初始化的终端乘子集合。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
void initializeFinalMultiplierCollection(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>& ocp,
    Scalar time,
    MultiplierCollection<Scalar, FinalStageConstraintLayout<ConstraintConfig>>&
        multiplierCollection) {
  ocp.finalEqualityLagrangian.initializeLagrangian(
      time, multiplierCollection.stateEq);
  ocp.finalInequalityLagrangian.initializeLagrangian(
      time, multiplierCollection.stateIneq);
}

/**
 * @brief 初始化中间时刻的乘子集合（等式与不等式拉格朗日）。
 * @param [in] ocp 最优控制问题。
 * @param [in] time 中间时刻。
 * @param [out] multiplierCollection 待初始化的中间乘子集合。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
void initializeIntermediateMultiplierCollection(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>& ocp,
    Scalar time,
    MultiplierCollection<Scalar,
                         IntermediateStageConstraintLayout<ConstraintConfig>>&
        multiplierCollection) {
  ocp.stateEqualityLagrangian.initializeLagrangian(
      time, multiplierCollection.stateEq);
  ocp.stateInequalityLagrangian.initializeLagrangian(
      time, multiplierCollection.stateIneq);
  ocp.equalityLagrangian.initializeLagrangian(
      time, multiplierCollection.stateInputEq);
  ocp.inequalityLagrangian.initializeLagrangian(
      time, multiplierCollection.stateInputIneq);
}

/**
 * @brief 根据缓存的对偶解初始化对偶解：若缓存非空则插值，否则用 ocp
 * 的拉格朗日初始化。
 * @param [in] ocp 最优控制问题。
 * @param [in] primalSolution 原始解（时间轨迹）。
 * @param [in] cachedDualSolution 缓存的对偶解（用于插值）。
 * @param [out] dualSolution 待初始化的对偶解。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
void initializeDualSolution(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>& ocp,
    const PrimalSolution<Scalar, Transcription>& primalSolution,
    const DualSolution<Scalar, typename Transcription::Horizon,
                       ConstraintConfig>& cachedDualSolution,
    DualSolution<Scalar, typename Transcription::Horizon, ConstraintConfig>&
        dualSolution) {
  constexpr std::size_t PredictLength = Transcription::PredictLength;
  using IntermediateMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           IntermediateStageConstraintLayout<ConstraintConfig>>;

  dualSolution.timeTrajectory = primalSolution.timeTrajectory_;

  if (!cachedDualSolution.empty()) {
    // final
    dualSolution.final = cachedDualSolution.final;

    // intermediates
    for (size_t i = 0; i < PredictLength; i++) {
      const Scalar& time = primalSolution.timeTrajectory_[i];
      IntermediateMultiplierCollection_t& multipliers =
          dualSolution.intermediates[i];
      multipliers = getIntermediateDualSolutionAtTime(cachedDualSolution, time);
    }
  } else {
    // final
    initializeFinalMultiplierCollection(
        ocp, primalSolution.timeTrajectory_.back(), dualSolution.final);

    // intermediates
    for (size_t i = 0; i < PredictLength; i++) {
      const Scalar& time = primalSolution.timeTrajectory_[i];
      IntermediateMultiplierCollection_t& multipliers =
          dualSolution.intermediates[i];
      initializeIntermediateMultiplierCollection(ocp, time, multipliers);
    }
  }
}

/**
 * @brief 根据当前状态-输入与 ocp 的拉格朗日更新规则原地更新对偶解，并同步更新
 * problemMetrics 的惩罚项。
 * @param [in] ocp 最优控制问题。
 * @param [in] primalSolution 原始解。
 * @param [in,out] problemMetrics 问题指标（其惩罚项将随对偶解更新）。
 * @param [out] dualSolution 待更新的对偶解（引用）。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
void updateDualSolution(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>& ocp,
    const PrimalSolution<Scalar, Transcription>& primalSolution,
    ProblemMetrics<Scalar, Transcription, ConstraintConfig>& problemMetrics,
    DualSolutionRef<Scalar, typename Transcription::Horizon, ConstraintConfig>
        dualSolution) {
  constexpr int XDim = Transcription::XDim;
  constexpr int UDim = Transcription::UDim;
  constexpr std::size_t PredictLength = Transcription::PredictLength;

  using FinalMetrics_t = Metrics<Scalar, typename Transcription::Dims,
                                 FinalStageConstraintLayout<ConstraintConfig>>;
  using IntermediateMetrics_t =
      Metrics<Scalar, typename Transcription::Dims,
              IntermediateStageConstraintLayout<ConstraintConfig>>;
  using FinalMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           FinalStageConstraintLayout<ConstraintConfig>>;
  using IntermediateMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           IntermediateStageConstraintLayout<ConstraintConfig>>;

  // final
  {
    const Scalar& time = primalSolution.timeTrajectory_.back();
    const Vector<Scalar, XDim>& state = primalSolution.stateTrajectory_.back();
    FinalMetrics_t& metrics = problemMetrics.final;
    FinalMultiplierCollection_t& multipliers = dualSolution.final;
    updateFinalMultiplierCollection(ocp, time, state, metrics, multipliers);
  }

  // intermediates
  // static_assert(dualSolution.intermediates.size() ==
  // primalSolution.timeTrajectory_.size());
  // static_assert(problemMetrics.intermediates.size() ==
  // primalSolution.timeTrajectory_.size());

  for (size_t i = 0; i < PredictLength; i++) {
    const Scalar& time = primalSolution.timeTrajectory_[i];
    const Vector<Scalar, XDim>& state = primalSolution.stateTrajectory_[i];
    const Vector<Scalar, UDim>& input = primalSolution.inputTrajectory_[i];
    IntermediateMetrics_t& metrics = problemMetrics.intermediates[i];
    IntermediateMultiplierCollection_t& multipliers =
        dualSolution.intermediates[i];

    updateIntermediateMultiplierCollection(ocp, time, state, input, metrics,
                                           multipliers);
  }
}

/**
 * @brief 原地更新终端乘子集合（等式与不等式），并同步更新 metrics 的惩罚项。
 * @param [in] ocp 最优控制问题。
 * @param [in] time 终端时间。
 * @param [in] state 终端状态。
 * @param [in,out] metrics 终端 Metrics（惩罚项将随乘子更新）。
 * @param [out] multipliers 待更新的终端乘子集合。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
void updateFinalMultiplierCollection(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>& ocp,
    Scalar time, const Vector<Scalar, Transcription::XDim>& state,
    Metrics<Scalar, typename Transcription::Dims,
            FinalStageConstraintLayout<ConstraintConfig>>& metrics,
    MultiplierCollection<Scalar, FinalStageConstraintLayout<ConstraintConfig>>&
        multipliers) {
  ocp.finalEqualityLagrangian.updateLagrangian(
      time, state, metrics.stateEqLagrangian, multipliers.stateEq);
  ocp.finalInequalityLagrangian.updateLagrangian(
      time, state, metrics.stateIneqLagrangian, multipliers.stateIneq);
}

/**
 * @brief 原地更新中间乘子集合（等式与不等式），并同步更新 metrics 的惩罚项。
 * @param [in] ocp 最优控制问题。
 * @param [in] time 中间时刻。
 * @param [in] state 中间状态。
 * @param [in] input 中间输入。
 * @param [in,out] metrics 中间 Metrics（惩罚项将随乘子更新）。
 * @param [out] multipliers 待更新的中间乘子集合。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
void updateIntermediateMultiplierCollection(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>& ocp,
    Scalar time, const Vector<Scalar, Transcription::XDim>& state,
    const Vector<Scalar, Transcription::UDim>& input,
    Metrics<Scalar, typename Transcription::Dims,
            IntermediateStageConstraintLayout<ConstraintConfig>>& metrics,
    MultiplierCollection<Scalar,
                         IntermediateStageConstraintLayout<ConstraintConfig>>&
        multipliers) {
  ocp.stateEqualityLagrangian.updateLagrangian(
      time, state, metrics.stateEqLagrangian, multipliers.stateEq);
  ocp.stateInequalityLagrangian.updateLagrangian(
      time, state, metrics.stateIneqLagrangian, multipliers.stateIneq);
  ocp.equalityLagrangian.updateLagrangian(time, state, input,
                                          metrics.stateInputEqLagrangian,
                                          multipliers.stateInputEq);
  ocp.inequalityLagrangian.updateLagrangian(time, state, input,
                                            metrics.stateInputIneqLagrangian,
                                            multipliers.stateInputIneq);
}
