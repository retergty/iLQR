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

#pragma once

#include "SearchStrategyBase.hpp"
#include "RolloutBase.hpp"
#include <functional>
#include "NumericTraits.hpp"

/**
 * Line search strategy: The class computes the nominal controller and the nominal trajectories as well the corresponding performance
 * indices. It line-searches on the feedforward parts of the controller and chooses the largest acceptable step-size.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
class LineSearchStrategy final : public SearchStrategyBase<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains,
                                                           FinalStateEqConstrains, FinalStateIneqConstrains>
{
public:
  using RolloutBase_t = RolloutBase<Scalar, XDimisions, UDimisions, PredictLength>;
  using OptimalControlProblem_t = OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using PerformanceIndex_t = PerformanceIndex<Scalar>;
  using DualSolution_t = DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength>;
  using SearchStrategySolution_t = SearchStrategySolution<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using SearchStrategySolutionRef_t = SearchStrategySolutionRef<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using StateVector_t = Vector<Scalar, XDimisions>;
  using InputVector_t = Vector<Scalar, UDimisions>;
  using LinearController_t = LinearController<Scalar, XDimisions, UDimisions, PredictLength>;
  using SearchStrategyBase_t = SearchStrategyBase<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains,
                                                  FinalStateEqConstrains, FinalStateIneqConstrains>;
  using ModelData_t = ModelData<Scalar, XDimisions, UDimisions>;
  /**
   * Constructor.
   *
   * @param [in] baseSettings: The basic settings for the search strategy algorithms.
   * @param [in] settings: The line search settings.
   * @param [in] threadPoolRef: A reference to the thread pool instance.
   * @param [in] rolloutRefStock: An array of references to the rollout.
   * @param [in] optimalControlProblemRef: An array of references to the optimal control problem.
   * @param [in] meritFunc: the merit function which gets the PerformanceIndex and returns the merit function value.
   */
  LineSearchStrategy(RolloutBase_t &rollout, OptimalControlProblem_t &optimalControlProblem,
                     std::function<Scalar(const PerformanceIndex_t &)> meritFunc) : SearchStrategyBase_t(), rollout_(rollout), optimalControlProblem_(optimalControlProblem), meritFunc_(meritFunc)
  {
    // Initialize controller
    workersSolution_.primalSolution.controllerPtr_.reset(new LinearController_t);
  }

  ~LineSearchStrategy() override = default;
  LineSearchStrategy(const LineSearchStrategy &) = delete;
  LineSearchStrategy &operator=(const LineSearchStrategy &) = delete;

  void reset() override {}

  bool run(const std::pair<Scalar, Scalar> &timePeriod, const StateVector_t &initState, const Scalar expectedCost,
           const LinearController_t &unoptimizedController, const DualSolution_t &dualSolution, SearchStrategySolutionRef_t &solutionRef) override
  {
    // initialize lineSearchModule inputs
    lineSearchInputRef_.timePeriodPtr = &timePeriod;
    lineSearchInputRef_.initStatePtr = &initState;
    lineSearchInputRef_.unoptimizedControllerPtr = &unoptimizedController;
    lineSearchInputRef_.dualSolutionPtr = &dualSolution;
    bestSolutionRef_ = &solutionRef;

    // perform a rollout with steplength zero.
    constexpr size_t taskId = 0;
    constexpr Scalar stepLength = 0.0;

    computeSolution(taskId, stepLength, workersSolution_);
    baselineMerit_ = workersSolution_.performanceIndex.merit;
    unoptimizedControllerUpdateIS_ = computeControllerUpdateIS(unoptimizedController);

    // record solution
    bestStepSize_ = stepLength;
    swap(*bestSolutionRef_, workersSolution_);

    lineSearchTask();

    return true;
  }

  bool checkConvergence(bool unreliableControllerIncrement,
                        const PerformanceIndex_t &previousPerformanceIndex,
                        const PerformanceIndex_t &currentPerformanceIndex) const override
  {
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

  void computeRiccatiModification(const ModelData_t &projectedModelData, Matrix<Scalar, XDimisions, XDimisions> &deltaQm) const override
  {
    // const auto &QmProjected = projectedModelData.cost.dfdxx;
    // const auto &PmProjected = projectedModelData.cost.dfdux;

    // Q_minus_PTRinvP
    // matrix_t Q_minus_PTRinvP = QmProjected;
    // Q_minus_PTRinvP.noalias() -= PmProjected.transpose() * PmProjected;

    // deltaQm
    // deltaQm = Q_minus_PTRinvP;
    deltaQm.setZero();
    shiftHessian(settings_.hessianCorrectionStrategy, deltaQm, settings_.hessianCorrectionMultiple);
    // deltaQm -= Q_minus_PTRinvP;
  }

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
  constexpr size_t maxNumOfSearches() const
  {
    size_t maxNumOfLineSearches = 0;
    if (numerics::almost_eq(settings_.minStepLength, settings_.maxStepLength))
    {
      maxNumOfLineSearches = 1;
    }
    else if (settings_.maxStepLength < settings_.minStepLength)
    {
      maxNumOfLineSearches = 0;
    }
    else
    {
      constexpr Scalar ratio = settings_.minStepLength / settings_.maxStepLength;
      // maxNumOfLineSearches =
      //     static_cast<size_t>(std::log(ratio + numeric_traits::limitEpsilon<Scalar>()) / std::log(settings_.contractionRate) + 1);
      // search it constexpr, find a maxNumOfLineSearches to satisfy settings_.maxStepLength * settings_.contractionRate^maxNumOfLineSearches < settings_.minStepLength
      Scalar result = 1;

      while (result >= ratio)
      {
        result *= settings_.contractionRate;
        maxNumOfLineSearches++;
      }
    }
    return maxNumOfLineSearches;
  }

  /** Computes the solution on a given stepLength  */
  void computeSolution(Scalar stepLength, SearchStrategySolution_t &solution)
  {
    // compute primal solution
    incrementController(stepLength, *lineSearchInputRef_.unoptimizedControllerPtr, getLinearController(solution.primalSolution));
    solution.avgTimeStep = rolloutTrajectory(rollout_, lineSearchInputRef_.timePeriodPtr->first, *lineSearchInputRef_.initStatePtr,
                                             lineSearchInputRef_.timePeriodPtr->second, solution.primalSolution);

    // adjust dual solution only if it is required
    const DualSolution_t *adjustedDualSolutionPtr = lineSearchInputRef_.dualSolutionPtr;

    // initialize dual solution
    initializeDualSolution(optimalControlProblem_, solution.primalSolution, *adjustedDualSolutionPtr, solution.dualSolution);

    // compute problem metrics
    computeRolloutMetrics(optimalControlProblem_, solution.primalSolution, solution.dualSolution, solution.problemMetrics);

    // compute performanceIndex
    solution.performanceIndex = computeRolloutPerformanceIndex(solution.primalSolution.timeTrajectory_, solution.problemMetrics);
    solution.performanceIndex.merit = meritFunc_(solution.performanceIndex);
  }

  /**
   * Defines line search task on a thread with various learning rates and choose the largest acceptable step-size.
   * The class computes the nominal controller and the nominal trajectories as well the corresponding performance indices.
   */
  void lineSearchTask()
  {
    Scalar stepLength = settings_.maxStepLength;

    // run workers
    for (int i = 0; i < alphaProcessed_.size(); ++i)
    {
      alphaProcessed_[i] = false;
    }

    for (int alphaExp = 0; alphaExp < maxNumOfSearches(); ++alphaExp)
    {
      /*
       * finish this thread's task since the learning rate is less than the minimum learning rate.
       * This means that the all the line search tasks are already processed or they are under
       * process in other threads.
       */

      computeSolution(stepLength, workersSolution_);

      // whether to accept the step or reject it
      bool terminateLinesearchTasks = false;
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
        swap(*bestSolutionRef_, workersSolution_);
        break;
      }
      alphaProcessed_[alphaExp] = true;

      stepLength *= settings_.contractionRate;
    }
  }

  constexpr static LineSearchSettings<Scalar> settings_{};

  DualSolution_t tempDualSolutions_;
  SearchStrategySolution_t workersSolution_;

  RolloutBase_t &rollout_;
  OptimalControlProblem_t &optimalControlProblem_;
  std::function<Scalar(const PerformanceIndex_t &)> meritFunc_;

  // input
  LineSearchInputRef lineSearchInputRef_;
  // output
  std::atomic<Scalar> bestStepSize_{0.0};
  SearchStrategySolutionRef_t *bestSolutionRef_;

  // convergence check
  Scalar baselineMerit_ = 0.0;                 // the merit of the rollout for zero learning rate
  Scalar unoptimizedControllerUpdateIS_ = 0.0; // integral of the squared (IS) norm of the controller update.

  // threading
  // std::atomic_size_t nextTaskId_{0};
  std::array<uint8_t, maxNumOfSearches()> alphaProcessed_;
};
