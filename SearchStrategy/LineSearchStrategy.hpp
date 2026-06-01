/**
 * @file LineSearchStrategy.hpp
 * @brief 线搜索策略：在控制器前馈上做 Armijo
 * 回溯，选取最大可接受步长并更新优化解。
 */
#pragma once
#include <algorithm>
#include <utility>

#include "Misc/Numerics.hpp"
#include "SearchStrategy/SearchStrategyBase.hpp"
#include "SearchStrategy/SearchStrategySettings.hpp"

template <typename Descriptor>
class iLQR;

/**
 * @brief 线搜索策略：用未优化控制器做 rollout，在步长上做 Armijo
 * 回溯，选取满足下降条件的最大步长并写回解。
 */
template <typename Descriptor>
class LineSearchStrategy final : public SearchStrategyBase<Descriptor> {
 public:
  using Types = iLQRTypes<Descriptor>;
  using Scalar = typename Types::Scalar;
  using LinearController_t = typename Types::LinearController_t;
  using StateVector_t = typename Types::StateVector_t;
  using InputVector_t = typename Types::InputVector_t;
  using PerformanceIndex_t = typename Types::PerformanceIndex_t;
  using ModelData_t = typename Types::ModelData_t;
  using DualSolution_t = typename Types::DualSolution_t;
  using SearchStrategySolution_t = SearchStrategySolution<Descriptor>;
  using SearchStrategySolutionRef_t = SearchStrategySolutionRef<Descriptor>;
  using HmMatrix_t = typename Types::HmMatrix_t;
  using SmMatrix_t = typename Types::SmMatrix_t;
  using iLQR_t = iLQR<Descriptor>;
  using RolloutBase_t = typename Types::RolloutBase_t;
  using OptimalControlProblem_t = typename Types::OptimalControlProblem_t;

  using SearchStrategyBase_t = SearchStrategyBase<Descriptor>;

  /** @brief 构造线搜索策略，绑定 iLQR 实例（用于 rollout、merit 等）。 */
  LineSearchStrategy(const LineSearchSettings<Scalar>& setting, iLQR_t& ilqr)
      : settings_(setting), ilqr_(ilqr) {}

  ~LineSearchStrategy() override = default;
  LineSearchStrategy(const LineSearchStrategy&) = delete;
  LineSearchStrategy& operator=(const LineSearchStrategy&) = delete;

  void reset() override {}

  /**
   * @brief 执行线搜索：从最大步长起尝试，选满足 Armijo
   * 的最大步长，将最优候选解写回 solutionRef。
   * @param [in] timePeriod 积分时间区间 (initTime, finalTime)。
   * @param [in] initState 初始状态。
   * @param [in] expectedCost 期望代价（当前未使用）。
   * @param [in] unoptimizedController 未优化控制器（含 deltaBias）。
   * @param [in] dualSolution 对偶解。
   * @param [in,out] solutionRef 输出解引用；primal/metrics/performanceIndex
   * 来自最佳候选解，dualSolution 复制自输入 dualSolution。
   * @return 当前实现恒返回 true。
   */
  bool run(const std::pair<Scalar, Scalar>& timePeriod,
           const StateVector_t& initState, const Scalar expectedCost,
           const LinearController_t& unoptimizedController,
           const DualSolution_t& dualSolution,
           SearchStrategySolutionRef_t& solutionRef) override {
    (void)expectedCost;
    const LineSearchInputRef inputRef{&timePeriod, &initState,
                                      &unoptimizedController, &dualSolution};

    SearchStrategySolution_t& baselineSolution = baselineSolutionCache_;
    SearchStrategySolution_t& candidateSolution = candidateSolutionCache_;

    initializeControllerStructure(unoptimizedController,
                                  baselineSolution.primalSolution.controller_);
    initializeControllerStructure(unoptimizedController,
                                  candidateSolution.primalSolution.controller_);

    // 以零步长真实 rollout 作为 Armijo 基准。
    computeSolution(inputRef, Scalar(0.0), baselineSolution);
    const Scalar baselineMerit = baselineSolution.performanceIndex.merit;
    const Scalar unoptimizedControllerUpdateIS =
        iLQR_t::computeControllerUpdateIS(unoptimizedController);

    const bool acceptedCandidate =
        lineSearchTask(inputRef, baselineMerit, unoptimizedControllerUpdateIS,
                       candidateSolution);
    const SearchStrategySolution_t& selectedSolution =
        acceptedCandidate ? candidateSolution : baselineSolution;

    solutionRef.primalSolution = selectedSolution.primalSolution;
    solutionRef.dualSolution = dualSolution;
    solutionRef.avgTimeStep = selectedSolution.avgTimeStep;
    solutionRef.performanceIndex = selectedSolution.performanceIndex;
    solutionRef.problemMetrics = selectedSolution.problemMetrics;
    return true;
  }

  /**
   * @brief 判断是否收敛：基于前后两次 performance index 的总代价相对变化与
   * minRelCost 比较。
   * @param [in] unreliableControllerIncrement 当前未使用。
   * @param [in] previousPerformanceIndex 上一迭代的 性能指标。
   * @param [in] currentPerformanceIndex 当前迭代的 性能指标。
   * @return 若总代价相对变化 <= minRelCost 则返回 true。
   */
  bool checkConvergence(
      bool unreliableControllerIncrement,
      const PerformanceIndex_t& previousPerformanceIndex,
      const PerformanceIndex_t& currentPerformanceIndex) const override {
    (void)unreliableControllerIncrement;
    // 循环中断变量。
    const Scalar currentTotalCost =
        currentPerformanceIndex.cost +
        currentPerformanceIndex.equalityLagrangian +
        currentPerformanceIndex.inequalityLagrangian;
    const Scalar previousTotalCost =
        previousPerformanceIndex.cost +
        previousPerformanceIndex.equalityLagrangian +
        previousPerformanceIndex.inequalityLagrangian;
    const Scalar absCostChange = std::abs(currentTotalCost - previousTotalCost);
    const Scalar costReference =
        std::max(settings_.costNormalizationBase, std::abs(previousTotalCost));
    const Scalar relCostChange = absCostChange / costReference;
    const bool isCostFunctionConverged = relCostChange <= settings_.minRelCost;
    const bool isOptimizationConverged = isCostFunctionConverged;

    return isOptimizationConverged;
  }

  /**
   * @brief 计算 Riccati 修正：当前实现将 deltaQm 置零并施加 Hessian
   * 修正（shift）。
   * @param [in] modelData 当前节点模型数据（当前未使用）。
   * @param [out] deltaQm 输出的 Riccati 修正矩阵（被设为 shift 后的零矩阵）。
   */
  void computeRiccatiModification(const ModelData_t& modelData,
                                  SmMatrix_t& deltaQm) const override {
    (void)modelData;
    (void)deltaQm;
    deltaQm.setZero();
    shiftHessian(settings_.hessianCorrectionStrategy, deltaQm,
                 settings_.hessianCorrectionMultiple);
  }

  /** @brief 对哈密顿量 Hessian 的额外修正；当前实现直接返回 Hm，不做修改。 */
  HmMatrix_t augmentHamiltonianHessian(const ModelData_t& modelData,
                                       const HmMatrix_t& Hm) const override {
    (void)modelData;
    return Hm;
  }

 private:
  struct LineSearchInputRef {
    const std::pair<Scalar, Scalar>* timePeriodPtr;
    const StateVector_t* initStatePtr;
    const LinearController_t* unoptimizedControllerPtr;
    const DualSolution_t* dualSolutionPtr;
  };

  static void initializeControllerStructure(const LinearController_t& source,
                                            LinearController_t& target) {
    target.timeStamp_ = source.timeStamp_;
    target.gainArray_ = source.gainArray_;
  }

  /** 线搜索迭代次数（if 语句顺序很重要）。 */
  constexpr static size_t maxNumOfSearches(
      const LineSearchSettings<Scalar> settings) {
    size_t maxNumOfLineSearches = 0;
    if (numerics::almost_eq(settings.minStepLength, settings.maxStepLength)) {
      maxNumOfLineSearches = 1;
    } else if (settings.maxStepLength < settings.minStepLength) {
      maxNumOfLineSearches = 0;
    } else {
      Scalar ratio = settings.minStepLength / settings.maxStepLength;
      Scalar result = 1;

      while (result >= ratio) {
        result *= settings.contractionRate;
        maxNumOfLineSearches++;
      }
    }
    return maxNumOfLineSearches;
  }

  /**
   * @brief 对给定步长计算候选解：更新控制器 bias、执行 rollout、计算 metrics 与
   * 性能指标。
   * @param [in] stepLength 线搜索步长。
   * @param [out] solution 输出的候选轨迹、metrics、performanceIndex。
   * dualSolution 只作为 metrics 计算输入，不存入候选解。
   */
  void computeSolution(const LineSearchInputRef& inputRef, Scalar stepLength,
                       SearchStrategySolution_t& solution) {
    // 计算原始解。
    iLQR_t::changeControllerStepLength(stepLength,
                                       *inputRef.unoptimizedControllerPtr,
                                       solution.primalSolution.controller_);
    solution.avgTimeStep = iLQR_t::rolloutTrajectory(
        ilqr_.rollout(), inputRef.timePeriodPtr->first, *inputRef.initStatePtr,
        inputRef.timePeriodPtr->second, solution.primalSolution);

    // 计算问题指标。
    iLQR_t::computeRolloutMetrics(
        ilqr_.optimalControlProblem(), ilqr_.targetTrajectory(),
        solution.primalSolution, *inputRef.dualSolutionPtr,
        solution.problemMetrics);

    // 计算 performanceIndex。
    solution.performanceIndex = iLQR_t::computeRolloutPerformanceIndex(
        solution.primalSolution.timeTrajectory_, solution.problemMetrics);
    solution.performanceIndex.merit =
        ilqr_.calculateRolloutMerit(solution.performanceIndex);
  }

  /** @brief 从最大步长开始按收缩率递减尝试，返回是否接受正步长候选。 */
  bool lineSearchTask(const LineSearchInputRef& inputRef, Scalar baselineMerit,
                      Scalar unoptimizedControllerUpdateIS,
                      SearchStrategySolution_t& candidateSolution) {
    Scalar stepLength = settings_.maxStepLength;
    const size_t MaxSearch = maxNumOfSearches(settings_);
    for (size_t alphaExp = 0; alphaExp < MaxSearch; ++alphaExp) {
      computeSolution(inputRef, stepLength, candidateSolution);

      /*
       * 基于 “Armijo backtracking” 步长选择策略：
       * 代价应优于基准代价，同时学习率应
       * 尽可能大。这等价于单核线搜索。
       * 搜索。
       */
      const bool armijoCondition =
          candidateSolution.performanceIndex.merit <
          (baselineMerit - settings_.armijoCoefficient * stepLength *
                               unoptimizedControllerUpdateIS);
      if (armijoCondition) {
        return true;
      }

      stepLength *= settings_.contractionRate;
    }
    return false;
  }

  LineSearchSettings<Scalar> settings_{};

  iLQR_t& ilqr_;

  SearchStrategySolution_t baselineSolutionCache_;
  SearchStrategySolution_t candidateSolutionCache_;
};
