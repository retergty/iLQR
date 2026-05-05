/******************************************************************************
Copyright (c) 2017, Farbod Farshidian. All rights reserved.

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
 * @file LineSearchStrategy.hpp
 * @brief 线搜索策略：在控制器前馈上做 Armijo 回溯，选取最大可接受步长并更新优化解。
 */
#pragma once

#include "SearchStrategyBase.hpp"
#include "RolloutBase.hpp"
#include <functional>
#include "NumericTraits.hpp"

template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
class iLQR;

/**
 * @brief 线搜索策略：用未优化控制器做 rollout，在步长上做 Armijo 回溯，选取满足下降条件的最大步长并写回解。
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
class LineSearchStrategy final : public SearchStrategyBase<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains,
                                                           FinalStateEqConstrains, FinalStateIneqConstrains>
{
public:
  using iLQR_t = iLQR<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using RolloutBase_t = RolloutBase<Scalar, XDimisions, UDimisions>;
  using OptimalControlProblem_t = OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using PerformanceIndex_t = PerformanceIndex<Scalar>;
  using DualSolution_t = DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength>;
  using SearchStrategySolution_t = SearchStrategySolution<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using SearchStrategySolutionRef_t = SearchStrategySolutionRef<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using StateVector_t = Vector<Scalar, XDimisions>;
  using InputVector_t = Vector<Scalar, UDimisions>;
  using LinearController_t = LinearController<Scalar, XDimisions, UDimisions, PredictLength + 1>;
  using SearchStrategyBase_t = SearchStrategyBase<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains,
                                                  FinalStateEqConstrains, FinalStateIneqConstrains>;
  using ModelData_t = ModelData<Scalar, XDimisions, UDimisions>;

  /** @brief 构造线搜索策略，绑定 iLQR 实例（用于 rollout、merit 等）。 */
  LineSearchStrategy(iLQR_t &ilqr) : ilqr_(ilqr)
  {
  }

  ~LineSearchStrategy() override = default;
  LineSearchStrategy(const LineSearchStrategy &) = delete;
  LineSearchStrategy &operator=(const LineSearchStrategy &) = delete;

  void reset() override {}

  /**
   * @brief 执行线搜索：从最大步长起尝试，选满足 Armijo 的最大步长，将最优轨迹与 dual 写回 solutionRef。
   * @param [in] timePeriod 积分时间区间 (initTime, finalTime)。
   * @param [in] initState 初始状态。
   * @param [in] expectedCost 期望代价（当前未使用）。
   * @param [in] unoptimizedController 未优化控制器（含 deltaBias）。
   * @param [in] dualSolution 对偶解。
   * @param [in,out] solutionRef 输出解引用（primal/dual/metrics/performanceIndex 被写入）。
   * @return 当前实现恒返回 true。
   */
  bool run(const std::pair<Scalar, Scalar> &timePeriod, const StateVector_t &initState, const Scalar expectedCost,
           const LinearController_t &unoptimizedController, const DualSolution_t &dualSolution, SearchStrategySolutionRef_t &solutionRef) override
  {
    (void)expectedCost;
    // initialize lineSearchModule inputs
    lineSearchInputRef_.timePeriodPtr = &timePeriod;
    lineSearchInputRef_.initStatePtr = &initState;
    lineSearchInputRef_.unoptimizedControllerPtr = &unoptimizedController;
    lineSearchInputRef_.dualSolutionPtr = &dualSolution;
    bestSolutionRef_ = &solutionRef;

    // perform a rollout with steplength zero.
    Scalar stepLength = 0.0;

    iLQR_t::incrementController(stepLength, *lineSearchInputRef_.unoptimizedControllerPtr, workersSolution_.primalSolution.controller_);

    computeSolution(stepLength, workersSolution_);
    baselineMerit_ = workersSolution_.performanceIndex.merit;
    unoptimizedControllerUpdateIS_ = iLQR_t::computeControllerUpdateIS(unoptimizedController);

    // record solution
    bestStepSize_ = stepLength;
    bestSolutionRef_->swap(workersSolution_);
    lineSearchTask();

    return true;
  }

  /**
   * @brief 判断是否收敛：基于前后两次 performance index 的总代价相对变化与 minRelCost 比较。
   * @param [in] unreliableControllerIncrement 当前未使用。
   * @param [in] previousPerformanceIndex 上一迭代的 performance index。
   * @param [in] currentPerformanceIndex 当前迭代的 performance index。
   * @return 若总代价相对变化 <= minRelCost 则返回 true。
   */
  bool checkConvergence(bool unreliableControllerIncrement,
                        const PerformanceIndex_t &previousPerformanceIndex,
                        const PerformanceIndex_t &currentPerformanceIndex) const override
  {
    (void)unreliableControllerIncrement;
    // loop break variables
    const Scalar currentTotalCost =
        currentPerformanceIndex.cost + currentPerformanceIndex.equalityLagrangian + currentPerformanceIndex.inequalityLagrangian;
    const Scalar previousTotalCost =
        previousPerformanceIndex.cost + previousPerformanceIndex.equalityLagrangian + previousPerformanceIndex.inequalityLagrangian;
    const Scalar relCost = std::abs(currentTotalCost - previousTotalCost);
    const bool isCostFunctionConverged = relCost <= this->baseSettings_.minRelCost;
    const bool isOptimizationConverged = isCostFunctionConverged;

    return isOptimizationConverged;
  }

  /**
   * @brief 计算 Riccati 修正：当前实现将 deltaQm 置零并施加 Hessian 修正（shift）。
   * @param [in] projectedModelData 投影后模型数据（当前未使用）。
   * @param [out] deltaQm 输出的 Riccati 修正矩阵（被设为 shift 后的零矩阵）。
   */
  void computeRiccatiModification(const ModelData_t &projectedModelData, Matrix<Scalar, XDimisions, XDimisions> &deltaQm) const override
  {
    // const auto &QmProjected = projectedModelData.cost.dfdxx;
    // const auto &PmProjected = projectedModelData.cost.dfdux;
    (void)projectedModelData;
    (void)deltaQm;
    // Q_minus_PTRinvP
    // matrix_t Q_minus_PTRinvP = QmProjected;
    // Q_minus_PTRinvP.noalias() -= PmProjected.transpose() * PmProjected;

    // deltaQm
    // deltaQm = Q_minus_PTRinvP;
    deltaQm.setZero();
    shiftHessian(settings_.hessianCorrectionStrategy, deltaQm, settings_.hessianCorrectionMultiple);
    // deltaQm -= Q_minus_PTRinvP;
  }

  /** @brief 对哈密顿量 Hessian 的额外修正；当前实现直接返回 Hm，不做修改。 */
  Matrix<Scalar, UDimisions, UDimisions> augmentHamiltonianHessian(const ModelData_t & /*modelData*/, const Matrix<Scalar, UDimisions, UDimisions> &Hm) const override { return Hm; }

private:
  struct LineSearchInputRef
  {
    const std::pair<Scalar, Scalar> *timePeriodPtr;
    const StateVector_t *initStatePtr;
    const LinearController_t *unoptimizedControllerPtr;
    const DualSolution_t *dualSolutionPtr;
  };

  /** number of line search iterations (the if statements order is important) */
  constexpr static size_t maxNumOfSearches(const LineSearchSettings<Scalar> settings)
  {
    size_t maxNumOfLineSearches = 0;
    if (numerics::almost_eq(settings.minStepLength, settings.maxStepLength))
    {
      maxNumOfLineSearches = 1;
    }
    else if (settings.maxStepLength < settings.minStepLength)
    {
      maxNumOfLineSearches = 0;
    }
    else
    {
      Scalar ratio = settings.minStepLength / settings.maxStepLength;
      // maxNumOfLineSearches =
      //     static_cast<size_t>(std::log(ratio + numeric_traits::limitEpsilon<Scalar>()) / std::log(settings_.contractionRate) + 1);
      // search it constexpr, find a maxNumOfLineSearches to satisfy settings_.maxStepLength * settings_.contractionRate^maxNumOfLineSearches < settings_.minStepLength
      Scalar result = 1;

      while (result >= ratio)
      {
        result *= settings.contractionRate;
        maxNumOfLineSearches++;
      }
    }
    return maxNumOfLineSearches;
  }

  /**
   * @brief 对给定步长计算解：更新控制器 bias、执行 rollout、计算 metrics 与 performance index。
   * @param [in] stepLength 线搜索步长。
   * @param [out] solution 输出的轨迹、dual、metrics、performanceIndex。
   */
  void computeSolution(Scalar stepLength, SearchStrategySolution_t &solution)
  {
    // compute primal solution
    iLQR_t::changeControllerStepLength(stepLength, *lineSearchInputRef_.unoptimizedControllerPtr, solution.primalSolution.controller_);
    solution.avgTimeStep = iLQR_t::rolloutTrajectory(ilqr_.rollout_, lineSearchInputRef_.timePeriodPtr->first, *lineSearchInputRef_.initStatePtr,
                                                     lineSearchInputRef_.timePeriodPtr->second, solution.primalSolution);

    // initialize dual solution
    // initializeDualSolution(ilqr_.optimalControlProblem_, solution.primalSolution, *adjustedDualSolutionPtr, solution.dualSolution);
    solution.dualSolution = *lineSearchInputRef_.dualSolutionPtr;

    // compute problem metrics
    iLQR_t::computeRolloutMetrics(ilqr_.optimalControlProblem_, solution.primalSolution, solution.dualSolution, solution.problemMetrics);

    // compute performanceIndex
    solution.performanceIndex = iLQR_t::computeRolloutPerformanceIndex(solution.primalSolution.timeTrajectory_, solution.problemMetrics);
    solution.performanceIndex.merit = ilqr_.calculateRolloutMerit(solution.performanceIndex);
  }

  /**
   * @brief 从最大步长开始按收缩率递减尝试，选取满足 Armijo 条件的最大步长并更新 bestSolutionRef_。
   */
  void lineSearchTask()
  {
    Scalar stepLength = settings_.maxStepLength;
    const size_t MaxSearch = maxNumOfSearches(settings_);
    for (size_t alphaExp = 0; alphaExp < MaxSearch; ++alphaExp)
    {
      /*
       * finish this thread's task since the learning rate is less than the minimum learning rate.
       * This means that the all the line search tasks are already processed or they are under
       * process in other threads.
       */

      computeSolution(stepLength, workersSolution_);

      /*
       * based on the "Armijo backtracking" step length selection policy:
       * cost should be better than the baseline cost but learning rate should
       * be as high as possible. This is equivalent to a single core line search.
       */
      const bool armijoCondition = workersSolution_.performanceIndex.merit <
                                   (baselineMerit_ - settings_.armijoCoefficient * stepLength * unoptimizedControllerUpdateIS_);
      if (armijoCondition && stepLength > bestStepSize_)
      { // save solution
        bestStepSize_ = stepLength;
        bestSolutionRef_->swap(workersSolution_);
        break;
      }

      stepLength *= settings_.contractionRate;
    }
  }

  constexpr static LineSearchSettings<Scalar> settings_{};

  iLQR_t &ilqr_;

  DualSolution_t tempDualSolutions_;
  SearchStrategySolution_t workersSolution_;

  // input
  LineSearchInputRef lineSearchInputRef_;
  // output
  std::atomic<Scalar> bestStepSize_{0.0};
  SearchStrategySolutionRef_t *bestSolutionRef_;

  // convergence check
  Scalar baselineMerit_ = 0.0;                 // the merit of the rollout for zero learning rate
  Scalar unoptimizedControllerUpdateIS_ = 0.0; // integral of the squared (IS) norm of the controller update.
};
