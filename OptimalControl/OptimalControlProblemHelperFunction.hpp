/**
 * @file OptimalControlProblemHelperFunction.hpp
 * @brief 最优控制问题辅助函数：初始化乘子、对偶解及近似等。
 */
#pragma once
#include "OptimalControlProblem.hpp"
#include "PrimalSolution.hpp"
#include "DualSolution.hpp"
#include "ProblemMetrics.hpp"
#include "Metrics.hpp"

/**
 * @brief 初始化终端时刻的乘子集合（等式与不等式拉格朗日）。
 * @param [in] ocp 最优控制问题。
 * @param [in] time 终端时间。
 * @param [out] multiplierCollection 待初始化的终端乘子集合。
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void initializeFinalMultiplierCollection(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                                         Scalar time, MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &multiplierCollection)
{
  ocp.finalEqualityLagrangian.initializeLagrangian(time, multiplierCollection.stateEq);
  ocp.finalInequalityLagrangian.initializeLagrangian(time, multiplierCollection.stateIneq);
}

/**
 * @brief 初始化中间时刻的乘子集合（等式与不等式拉格朗日）。
 * @param [in] ocp 最优控制问题。
 * @param [in] time 中间时刻。
 * @param [out] multiplierCollection 待初始化的中间乘子集合。
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void initializeIntermediateMultiplierCollection(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                                                Scalar time, MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multiplierCollection)
{
  ocp.stateEqualityLagrangian.initializeLagrangian(time, multiplierCollection.stateEq);
  ocp.stateInequalityLagrangian.initializeLagrangian(time, multiplierCollection.stateIneq);
  ocp.equalityLagrangian.initializeLagrangian(time, multiplierCollection.stateInputEq);
  ocp.inequalityLagrangian.initializeLagrangian(time, multiplierCollection.stateInputIneq);
}

/**
 * @brief 根据缓存的对偶解初始化对偶解：若缓存非空则插值，否则用 ocp 的拉格朗日初始化。
 * @param [in] ocp 最优控制问题。
 * @param [in] primalSolution 原始解（时间轨迹）。
 * @param [in] cachedDualSolution 缓存的对偶解（用于插值）。
 * @param [out] dualSolution 待初始化的对偶解。
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void initializeDualSolution(
    const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
    const PrimalSolution<Scalar, XDimisions, UDimisions, PredictLength> &primalSolution,
    const DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength> &cachedDualSolution,
    DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength> &dualSolution)
{
  dualSolution.timeTrajectory = primalSolution.timeTrajectory_;

  if (!cachedDualSolution.empty())
  {
    // final
    dualSolution.final = cachedDualSolution.final;
  }
  else
  {
    initializeFinalMultiplierCollection(ocp, primalSolution.timeTrajectory_.back(), dualSolution.final);
  }

  if (!cachedDualSolution.empty())
  {
    // intermediates
    for (size_t i = 0; i < PredictLength; i++)
    {
      const Scalar &time = primalSolution.timeTrajectory_[i];
      MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multipliers = dualSolution.intermediates[i];
      multipliers = getIntermediateDualSolutionAtTime(cachedDualSolution, time);
    }
  }
  else
  {
    // intermediates
    for (size_t i = 0; i < PredictLength; i++)
    {
      const Scalar &time = primalSolution.timeTrajectory_[i];
      MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multipliers = dualSolution.intermediates[i];
      initializeIntermediateMultiplierCollection(ocp, time, multipliers);
    }
  }
}

/**
 * @brief 根据当前状态-输入与 ocp 的拉格朗日更新规则原地更新对偶解，并同步更新 problemMetrics 的惩罚项。
 * @param [in] ocp 最优控制问题。
 * @param [in] primalSolution 原始解。
 * @param [in,out] problemMetrics 问题指标（其惩罚项将随对偶解更新）。
 * @param [out] dualSolution 待更新的对偶解（引用）。
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void updateDualSolution(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                        const PrimalSolution<Scalar, XDimisions, UDimisions, PredictLength> &primalSolution,
                        ProblemMetrics<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &problemMetrics,
                        DualSolutionRef<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength> dualSolution)
{
  // final
  {
    const Scalar &time = primalSolution.timeTrajectory_.back();
    const Vector<Scalar, XDimisions> &state = primalSolution.stateTrajectory_.back();
    Metrics<Scalar, XDimisions, UDimisions, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &metrics = problemMetrics.final;
    MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &multipliers = dualSolution.final;
    updateFinalMultiplierCollection(ocp, time, state, metrics, multipliers);
  }

  // intermediates
  // static_assert(dualSolution.intermediates.size() == primalSolution.timeTrajectory_.size());
  // static_assert(problemMetrics.intermediates.size() == primalSolution.timeTrajectory_.size());

  for (size_t i = 0; i < PredictLength; i++)
  {
    const Scalar &time = primalSolution.timeTrajectory_[i];
    const Vector<Scalar, XDimisions> &state = primalSolution.stateTrajectory_[i];
    const Vector<Scalar, UDimisions> &input = primalSolution.inputTrajectory_[i];
    Metrics<Scalar, XDimisions, UDimisions, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &metrics = problemMetrics.intermediates[i];
    MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multipliers = dualSolution.intermediates[i];

    updateIntermediateMultiplierCollection(ocp, time, state, input, metrics, multipliers);
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
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void updateFinalMultiplierCollection(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                                     Scalar time, const Vector<Scalar, XDimisions> &state,
                                     Metrics<Scalar, XDimisions, UDimisions, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &metrics,
                                     MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &multipliers)
{
  ocp.finalEqualityLagrangian.updateLagrangian(time, state, metrics.stateEqLagrangian, multipliers.stateEq);
  ocp.finalInequalityLagrangian.updateLagrangian(time, state, metrics.stateIneqLagrangian, multipliers.stateIneq);
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
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void updateIntermediateMultiplierCollection(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                                            Scalar time, const Vector<Scalar, XDimisions> &state, const Vector<Scalar, UDimisions> &input,
                                            Metrics<Scalar, XDimisions, UDimisions, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &metrics,
                                            MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multipliers)
{
  ocp.stateEqualityLagrangian.updateLagrangian(time, state, metrics.stateEqLagrangian, multipliers.stateEq);
  ocp.stateInequalityLagrangian.updateLagrangian(time, state, metrics.stateIneqLagrangian, multipliers.stateIneq);
  ocp.equalityLagrangian.updateLagrangian(time, state, input, metrics.stateInputEqLagrangian, multipliers.stateInputEq);
  ocp.inequalityLagrangian.updateLagrangian(time, state, input, metrics.stateInputIneqLagrangian, multipliers.stateInputIneq);
}
