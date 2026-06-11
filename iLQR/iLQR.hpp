/**
 * @file iLQR.hpp
 * @brief 离散时间 iLQR 求解器：在名义轨迹上做 LQ 近似并迭代求解
 * Riccati，支持线搜索。
 */
#pragma once
#include <array>
#include <cassert>

#include "Integration/TrapezoidalIntegration.hpp"
#include "OptimalControl/OptimalControlProblemHelperFunction.hpp"
#include "OptimalControlData/PerformanceIndex.hpp"
#include "SearchStrategy/LineSearchStrategy.hpp"
#include "SearchStrategy/SearchStrategyBase.hpp"
#include "iLQR/DDPSetting.hpp"
#include "iLQR/HessianCorrection.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"
#include "iLQR/iLQRTypes.hpp"

/**
 * @brief 迭代线性二次调节器（iLQR）：基于名义轨迹的 LQ 近似与离散时间 Riccati
 * 反向递推的 DDP 求解器。
 * @tparam Descriptor
 * 编译期描述符，集中定义标量类型、维度、时域、动力学模式和约束配置。
 */
template <typename Descriptor>
class iLQR {
 public:
  using Types = iLQRTypes<Descriptor>;
  using Scalar = typename Types::Scalar;

  static constexpr int XDim = Types::XDim;
  static constexpr int UDim = Types::UDim;
  static constexpr std::size_t PredictLength = Types::PredictLength;

  using OptimalControlProblem_t = typename Types::OptimalControlProblem_t;

  using StateVector_t = typename Types::StateVector_t;
  using InputVector_t = typename Types::InputVector_t;
  using LvVector_t = typename Types::LvVector_t;
  using KmMatrix_t = typename Types::KmMatrix_t;
  using SmMatrix_t = typename Types::SmMatrix_t;
  using SvVector_t = typename Types::SvVector_t;
  using GmMatrix_t = typename Types::GmMatrix_t;
  using HmMatrix_t = typename Types::HmMatrix_t;
  using GvVector_t = typename Types::GvVector_t;

  using ModelData_t = typename Types::ModelData_t;
  using IntermediateMultiplierCollection_t =
      typename Types::IntermediateMultiplierCollection_t;
  using FinalMultiplierCollection_t =
      typename Types::FinalMultiplierCollection_t;
  using IntermediateMetrics_t = typename Types::IntermediateMetrics_t;
  using FinalMetrics_t = typename Types::FinalMetrics_t;
  using RolloutBase_t = typename Types::RolloutBase_t;
  using InitializerRollout_t = typename Types::InitializerRollout_t;
  using Initializer_t = typename Types::Initializer_t;

  using Rollout_t = typename Types::Rollout_t;
  using RolloutTrajectoryPointer_t = typename Types::RolloutTrajectoryPointer_t;

  using SearchStrategySolution_t = typename Types::SearchStrategySolution_t;
  using SearchStrategySolutionRef_t =
      typename Types::SearchStrategySolutionRef_t;

  using LinearController_t = typename Types::LinearController_t;
  using ControllerGainTrajectory_t = std::array<KmMatrix_t, PredictLength + 1>;
  using ControllerDeltaBiasTrajectory_t =
      std::array<LvVector_t, PredictLength + 1>;
  using SearchStrategyBase_t = typename Types::SearchStrategyBase_t;
  using LineSearchStrategy_t = LineSearchStrategy<Descriptor>;
  using RiccatiModification_t = typename Types::RiccatiModification_t;

  using TimeTrajectory_t = typename Types::TimeTrajectory_t;
  using StateTrajectory_t = typename Types::StateTrajectory_t;
  using InputTrajectory_t = typename Types::InputTrajectory_t;
  using TargetTrajectories_t = typename Types::TargetTrajectories_t;

  using IntermediateMultiplierTrajectory_t =
      typename Types::IntermediateMultiplierTrajectory_t;
  using ModelDataTrajectory_t = typename Types::ModelDataTrajectory_t;

  using PrimalSolution_t = typename Types::PrimalSolution_t;
  using DualSolution_t = typename Types::DualSolution_t;
  using DualSolutionRef_t = typename Types::DualSolutionRef_t;
  using LinearQuadraticApproximator_t =
      typename Types::LinearQuadraticApproximator_t;
  using PrimalDataContainer_t = typename Types::PrimalDataContainer_t;
  using DualDataContainer_t = typename Types::DualDataContainer_t;
  using ProblemMetrics_t = typename Types::ProblemMetrics_t;
  using PerformanceIndex_t = typename Types::PerformanceIndex_t;
  using PerformanceIndexEvaluator_t =
      typename Types::PerformanceIndexEvaluator_t;
  using IntermediatePerformanceIndexTrajectory_t =
      typename Types::IntermediatePerformanceIndexTrajectory_t;

  using Transcription_t = typename Types::Transcription_t;

  using ValueFunctionQuadraticApproximation_t =
      typename Types::ValueFunctionQuadraticApproximation_t;
  using ValueFunctionTrajectory_t = typename Types::ValueFunctionTrajectory_t;
  using DiscreteTimeRiccatiEquations_t =
      typename Types::DiscreteTimeRiccatiEquations_t;

  /**
   * @brief 构造 iLQR 求解器，绑定外部最优控制问题与初始化器。
   * @param [in] optimalControlProblem 最优控制问题定义，生命周期需长于求解器。
   * @param [in] initializer 轨迹初始化器，用于填补无控制器的时间段。
   */
  iLQR(const DDPSettings<Scalar>& ddp_setting,
       OptimalControlProblem_t& optimalControlProblem,
       Initializer_t* initializer)
      : ddpSettings_(ddp_setting),
        optimalControlProblem_(optimalControlProblem),
        rollout_(optimalControlProblem_.dynamicsPtr, ddpSettings_.timeStep),
        initializerRollout_(*initializer, ddpSettings_.timeStep),
        lineSearchStrategy_(makeLineSearchSettings(ddp_setting), *this) {
    assert(optimalControlProblem_.dynamicsPtr != nullptr);
    reset();
  };

  /**
   * @brief 重置求解器缓存，清除 warm start 解与运行统计。
   *
   * 不修改最优控制问题、参考轨迹、rollout 配置和临时工作区；下一次 run()
   * 将从 initializer 冷启动。
   */
  void reset() {
    initTime_ = Scalar(0.0);
    finalTime_ = Scalar(0.0);
    lastFinalTime_ = Scalar(0.0);

    optimizedPrimalSolution_.clear();
    optimizedDualSolution_.clear();

    totalNumIterations_ = 0;
    totalNumRuns_ = 0;

    lineSearchStrategy_.reset();
  }

  /**
   * 求解器主流程：针对给定初始
   * 状态、初始时间和终止时间运行优化器。
   *
   * @param [in] initTime 初始时间。
   * @param [in] initState 初始状态。
   */
  void run(Scalar initTime, const StateVector_t& initState) {
    // 初始化参数。
    initTime_ = initTime;
    initState_ = initState;
    lastFinalTime_ = finalTime_;
    finalTime_ = initTime + ddpSettings_.timeStep * (PredictLength);
    const size_t initIteration = totalNumIterations_;
    totalNumRuns_++;
    // optimized --> nominal：基于优化解初始化名义原始解和对偶解。
    // 基于优化解。
    bool initialSolutionExists =
        initializePrimalSolution();  // 如果 rollout 不完全来自
                                     // Initializer，则为 true。
    initializeDualSolutionAndMetrics(initialSolutionExists);

    performanceIndexLast_ = performanceIndex_;
    // 主循环收敛变量。
    bool isConverged = false;

    // DDP 主循环。
    while (true) {
      // nominal --> nominal：围绕名义轨迹构造 LQ 问题。
      // 轨迹。
      approximateOptimalControlProblem();

      // nominal --> nominal：求解 LQ 问题，直接写入待线搜索控制器的
      // feedback gain 与 feedforward delta bias。
      solveSequentialRiccatiEquations(
          nominalPrimalData_.modelDataFinalTime.cost,
          unoptimizedController_.gainArray_,
          unoptimizedController_.deltaBiasArray_);

      // 计算控制器并将结果存入 unoptimizedController_。
      calculateController();

      // Riccati 解计算出的期望代价/merit 不
      // 可靠。
      const auto lqModelExpectedCost =
          initialSolutionExists
              ? nominalDualData_.valueFunctionTrajectory.front().f
              : performanceIndex_.merit;

      // nominal --> optimized：基于当前 LQ 解更新
      // 优化后的原始解和对偶解。
      takePrimalDualStep(lqModelExpectedCost);

      // 迭代信息。
      ++totalNumIterations_;

      // 检查收敛。
      isConverged = lineSearchStrategy_.checkConvergence(
          !initialSolutionExists, performanceIndexLast_, performanceIndex_);
      performanceIndexLast_ = performanceIndex_;
      initialSolutionExists = true;

      if (isConverged || (totalNumIterations_ - initIteration) ==
                             ddpSettings_.maxNumIterations) {
        break;
      } else {
        // optimized --> nominal：将优化解作为下一次迭代的名义解。
        // 下一次迭代。
        optimizedDualSolution_.swap(nominalDualData_.dualSolution);
        optimizedPrimalSolution_.swap(nominalPrimalData_.primalSolution);
        optimizedProblemMetrics_.swap(nominalPrimalData_.problemMetrics);
      }
    }  // while 循环结束。
  }

  /**
   * @brief 设置参考轨迹（时间、状态、输入），用于跟踪型代价与初始化。
   * @param [in] timeTrajectory 时间序列，长度为 PredictLength+1。
   * @param [in] stateTrajectory 参考状态轨迹。
   * @param [in] inputTrajectory 参考输入轨迹。
   */
  void setDesireTrajectory(const TimeTrajectory_t& timeTrajectory,
                           const StateTrajectory_t& stateTrajectory,
                           const InputTrajectory_t& inputTrajectory) {
    targetTrajectory_.setTrajectory(timeTrajectory, stateTrajectory,
                                    inputTrajectory);
  }

  /**
   * @brief 设置下一次 run() 使用的固定时间步长。
   *
   * 同步更新 DDP 配置、主 rollout 和 initializer rollout 的步长。该接口保留
   * 已缓存的优化解，使下一次 run() 可以将旧控制器作为 warm start，并在新的
   * 时间网格上重新 rollout 与离散化。
   *
   * @param [in] dt 新时间步长，必须为正数。
   */
  void setTimeStep(const Scalar dt) {
    assert(dt > Scalar(0));
    ddpSettings_.timeStep = dt;
    rollout_.settings().timeStep = dt;
    initializerRollout_.settings().timeStep = dt;
  }

  const PrimalSolution_t& primalSolution() const {
    return optimizedPrimalSolution_;
  }
  const PerformanceIndex_t& performanceIndex() const {
    return performanceIndex_;
  }
  const DDPSettings<Scalar>& ddpSettings() const { return ddpSettings_; }
  size_t totalNumIterations() const { return totalNumIterations_; }
  size_t totalNumRuns() const { return totalNumRuns_; }
  Scalar averageNumIterations() const {
    return totalNumRuns_ == 0 ? Scalar(0)
                              : static_cast<Scalar>(totalNumIterations_) /
                                    static_cast<Scalar>(totalNumRuns_);
  }
  TargetTrajectories_t& targetTrajectory() { return targetTrajectory_; }
  const TargetTrajectories_t& targetTrajectory() const {
    return targetTrajectory_;
  }
  OptimalControlProblem_t& optimalControlProblem() {
    return optimalControlProblem_;
  }
  const OptimalControlProblem_t& optimalControlProblem() const {
    return optimalControlProblem_;
  }
  Rollout_t& rollout() { return rollout_; }
  const Rollout_t& rollout() const { return rollout_; }

  ValueFunctionQuadraticApproximation_t getValueFunction(
      const Scalar time, const StateVector_t& state) const {
    return getValueFunctionImpl(time, state, nominalPrimalData_.primalSolution,
                                nominalDualData_.valueFunctionTrajectory);
  }

  ValueFunctionQuadraticApproximation_t getValueFunction(
      size_t timeIndex, const StateVector_t& state) const {
    return getValueFunctionImpl(timeIndex, state,
                                nominalPrimalData_.primalSolution,
                                nominalDualData_.valueFunctionTrajectory);
  }

  ValueFunctionQuadraticApproximation_t getQFunction(
      size_t timeIndex, const StateVector_t& state,
      const InputVector_t& input) const {
    assert(timeIndex < PredictLength);

    const ModelData_t& modelData =
        nominalPrimalData_.modelDataTrajectory[timeIndex];
    const ValueFunctionQuadraticApproximation_t& nextValueFunction =
        nominalDualData_.valueFunctionTrajectory[timeIndex + 1];

    ValueFunctionQuadraticApproximation_t qFunction = modelData.cost;

    // Q_k(x,u) = l_k(x,u) + V_{k+1}(A dx + B du + h)
    const SmMatrix_t VxxA = nextValueFunction.dfdxx * modelData.dynamics.dfdx;
    const StateVector_t VxWithBias =
        nextValueFunction.dfdx + nextValueFunction.dfdxx * modelData.dynamics.f;

    qFunction.f +=
        nextValueFunction.f + nextValueFunction.dfdx.dot(modelData.dynamics.f) +
        Scalar(0.5) * modelData.dynamics.f.dot(nextValueFunction.dfdxx *
                                               modelData.dynamics.f);
    qFunction.dfdx += modelData.dynamics.dfdx.transpose() * VxWithBias;
    qFunction.dfdu += modelData.dynamics.dfdu.transpose() * VxWithBias;
    qFunction.dfdxx += modelData.dynamics.dfdx.transpose() * VxxA;
    qFunction.dfdux += modelData.dynamics.dfdu.transpose() * VxxA;
    qFunction.dfduu += modelData.dynamics.dfdu.transpose() *
                       nextValueFunction.dfdxx * modelData.dynamics.dfdu;

    const StateVector_t deltaX =
        state - nominalPrimalData_.primalSolution.stateTrajectory_[timeIndex];
    const InputVector_t deltaU =
        input - nominalPrimalData_.primalSolution.inputTrajectory_[timeIndex];
    const StateVector_t QxxDeltaX = qFunction.dfdxx * deltaX;
    const InputVector_t QuxDeltaX = qFunction.dfdux * deltaX;
    const InputVector_t QuuDeltaU = qFunction.dfduu * deltaU;

    qFunction.f +=
        deltaX.dot(Scalar(0.5) * QxxDeltaX + qFunction.dfdx) +
        deltaU.dot(QuxDeltaX + Scalar(0.5) * QuuDeltaU + qFunction.dfdu);
    qFunction.dfdx += QxxDeltaX + qFunction.dfdux.transpose() * deltaU;
    qFunction.dfdu += QuxDeltaX + QuuDeltaU;

    return qFunction;
  }

  /**
   * 根据性能指标计算 merit 函数。
   *
   * @param [in] performanceIndex 性能指标（含 cost、等式/不等式拉格朗日等）。
   * @return 用于线搜索与收敛判据的 merit 值。
   */
  static Scalar calculateRolloutMerit(
      const PerformanceIndex_t& performanceIndex) {
    // 代价
    Scalar merit = performanceIndex.cost;
    // 状态/状态-输入等式拉格朗日项。
    merit += performanceIndex.equalityLagrangian;
    // 状态/状态-输入不等式拉格朗日项。
    merit += performanceIndex.inequalityLagrangian;

    return merit;
  }

  /**
   * @brief 计算 rollout 轨迹中每个节点的代价、软约束和约束指标。
   *
   * @param [in] problem 最优控制问题。
   * @param [in] targetTrajectory 参考轨迹。
   * @param [in] primalSolution 原始解（轨迹）。
   * @param [in] dualSolution 对偶解。
   * @param [out] problemMetrics 各时刻代价、软约束与约束值。
   */
  static void computeRolloutMetrics(
      OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectory,
      const PrimalSolution_t& primalSolution,
      const DualSolution_t& dualSolution, ProblemMetrics_t& problemMetrics) {
    const TimeTrajectory_t& tTrajectory = primalSolution.timeTrajectory_;
    const StateTrajectory_t& xTrajectory = primalSolution.stateTrajectory_;
    const InputTrajectory_t& uTrajectory = primalSolution.inputTrajectory_;

    for (size_t k = 0; k < PredictLength; k++) {
      // 中间时刻代价和约束。
      problemMetrics.intermediates[k] =
          LinearQuadraticApproximator_t::computeIntermediateMetrics(
              problem, targetTrajectory, tTrajectory[k], xTrajectory[k],
              uTrajectory[k], dualSolution.intermediates[k]);
    }

    // 终端时刻代价和约束。
    problemMetrics.final = LinearQuadraticApproximator_t::computeFinalMetrics(
        problem, targetTrajectory, tTrajectory.back(), xTrajectory.back(),
        dualSolution.final);
  }

  /**
   * @brief 使用给定控制器前向展开系统轨迹。
   *
   * 连续动力学模式下，rollout 会积分系统流映射；离散动力学模式下，
   * rollout 会逐节点调用离散状态转移。
   *
   * @param [in] rollout 前向展开使用的 rollout 对象。
   * @param [in] initTime 初始时间。
   * @param [in] initState 初始状态。
   * @param [in] finalTime 终止时间。
   * @param [in,out] primalSolution 结果写入其时间/状态/输入轨迹；需已设置
   * controller 供 rollout 使用。
   */
  static int rolloutTrajectory(RolloutBase_t& rollout, Scalar initTime,
                               const StateVector_t& initState, Scalar finalTime,
                               PrimalSolution_t& primalSolution) {
    RolloutTrajectoryPointer_t rolloutTrajectortPtr(
        primalSolution.timeTrajectory_.data(),
        primalSolution.stateTrajectory_.data(),
        primalSolution.inputTrajectory_.data(), PredictLength + 1);
    return rollout.run(initTime, initState, finalTime,
                       &primalSolution.controller_, rolloutTrajectortPtr);
  }

  /**
   * 计算给定 ProblemMetrics 对应的 PerformanceIndex。
   *
   * @param [in] timeTrajectory rollout 的时间戳序列。
   * @param [in] problemMetrics 各时刻代价与约束指标。
   * @return 整条轨迹的 PerformanceIndex。
   */
  static PerformanceIndex_t computeRolloutPerformanceIndex(
      const TimeTrajectory_t& timeTrajectory,
      const ProblemMetrics_t& problemMetrics) {
    // 终端
    PerformanceIndex_t finalperformanceIndex =
        toPerformanceIndex(problemMetrics.final);
    // 中间节点
    IntermediatePerformanceIndexTrajectory_t performanceIndexTrajectory;

    for (size_t i = 0; i < performanceIndexTrajectory.size(); ++i) {
      performanceIndexTrajectory[i] =
          toPerformanceIndex(problemMetrics.intermediates[i]);
    }

    return PerformanceIndexEvaluator_t::evaluate(
        timeTrajectory, performanceIndexTrajectory, finalperformanceIndex);
  }

  /**
   * @brief 仅更新控制器的 bias：biasArray = unoptimizedController.biasArray +
   * stepLength * deltaBiasArray。
   * @param [in] stepLength 步长。
   * @param [in] unoptimizedController 源控制器。
   * @param [out] controller 目标控制器（仅 biasArray_ 被写入）。
   */
  static void changeControllerStepLength(
      Scalar stepLength, const LinearController_t& unoptimizedController,
      LinearController_t& controller) {
    for (size_t k = 0; k < unoptimizedController.size(); k++) {
      controller.biasArray_[k] =
          unoptimizedController.biasArray_[k] +
          stepLength * unoptimizedController.deltaBiasArray_[k];
    }
  }

  /**
   * 计算控制器更新平方范数的积分（IS）。
   *
   * @param [in] controller 控制器（含 deltaBiasArray）。
   * @return 控制器更新量的平方范数沿时间的积分（梯形积分）。
   */
  static Scalar computeControllerUpdateIS(
      const LinearController_t& controller) {
    std::array<Scalar, LinearController_t::size()> biasArraySquaredNorm;

    for (size_t i = 0; i < controller.size(); ++i) {
      biasArraySquaredNorm[i] = controller.deltaBiasArray_[i].squaredNorm();
    }
    // 使用梯形近似方法积分。
    return trapezoidalIntegration(controller.timeStamp_, biasArraySquaredNorm,
                                  Scalar(0.0));
  }

 private:
  static LineSearchSettings<Scalar> makeLineSearchSettings(
      const DDPSettings<Scalar>& ddpSettings) {
    LineSearchSettings<Scalar> settings = ddpSettings.lineSearch;
    settings.minRelCost = ddpSettings.minRelCost;
    return settings;
  }

  /**
   * @brief 基于上一轮优化解初始化当前名义原始解。
   *
   * @return 如果 rollout 不完全来自 Initializer，则返回 true。
   */
  bool initializePrimalSolution() {
    // 尝试使用控制器初始化。
    int numSteps = 0;
    bool ret = false;
    if (lastFinalTime_ > initTime_) {
      numSteps = rolloutInitialController(optimizedPrimalSolution_,
                                          nominalPrimalData_.primalSolution);
      ret = true;
    }
    // 距离上次运行已过去太久。
    rolloutInitializer(nominalPrimalData_.primalSolution, numSteps);
    return ret;
  }

  /**
   * @brief 使用 inputPrimalSolution 中的控制器前向展开系统轨迹。
   *
   * 通常会从当前初始状态展开到 finalTime_。如果 inputPrimalSolution
   * 的控制器没有覆盖完整区间，则只展开到该控制器可用的终止时间。
   *
   * @param [in] inputPrimalSolution 其控制器将用于前向 rollout。
   * @param [out] outputPrimalSolution 得到的原始解（轨迹与控制器）。
   * @return 被覆盖的步数。
   */
  int rolloutInitialController(PrimalSolution_t& inputPrimalSolution,
                               PrimalSolution_t& outputPrimalSolution) {
    outputPrimalSolution.controller_ = inputPrimalSolution.controller_;
    const Scalar rolloutFinalTime = std::min(lastFinalTime_, finalTime_);
    const int writtenPoints =
        rolloutTrajectory(rollout_, initTime_, initState_, rolloutFinalTime,
                          outputPrimalSolution);

    const int numSteps = writtenPoints - 1;
    assert(writtenPoints > 0);
    assert(numSteps <= static_cast<int>(PredictLength));

    return numSteps;
  }

  /**
   * 该过程会检查 primalSolution 的内容；如果其终止时间
   * 小于当前求解器的 finalTime_，则将其与
   * Initializer 的结果拼接。
   */
  void rolloutInitializer(PrimalSolution_t& primalSolution, int numSteps) {
    if (numSteps >= static_cast<int>(PredictLength)) {
      return;
    }

    RolloutTrajectoryPointer_t rolloutTrajectoryPtr(
        primalSolution.timeTrajectory_.data() + numSteps,
        primalSolution.stateTrajectory_.data() + numSteps,
        primalSolution.inputTrajectory_.data() + numSteps,
        PredictLength - numSteps + 1);
    if (numSteps == 0) {
      initializerRollout_.run(initTime_, initState_, finalTime_, nullptr,
                              rolloutTrajectoryPtr);
    } else {
      const Scalar initTime = primalSolution.timeTrajectory_[numSteps];
      const StateVector_t& initState =
          primalSolution.stateTrajectory_[numSteps];

      initializerRollout_.run(initTime, initState, finalTime_, nullptr,
                              rolloutTrajectoryPtr);
    }
  }

  /**
   * 基于优化解和名义原始解初始化名义对偶解，
   * 同时更新 ProblemMetrics。
   */
  void initializeDualSolutionAndMetrics(bool useCachedDualSolution) {
    // 初始化对偶解。
    initializeDualSolution(
        optimalControlProblem_, nominalPrimalData_.primalSolution,
        optimizedDualSolution_, nominalDualData_.dualSolution,
        useCachedDualSolution);

    computeRolloutMetrics(optimalControlProblem_, targetTrajectory_,
                          nominalPrimalData_.primalSolution,
                          nominalDualData_.dualSolution,
                          nominalPrimalData_.problemMetrics);

    // 计算 rollout merit。
    performanceIndex_ = computeRolloutPerformanceIndex(
        nominalPrimalData_.primalSolution.timeTrajectory_,
        nominalPrimalData_.problemMetrics);
    performanceIndex_.merit = calculateRolloutMerit(performanceIndex_);
  }

  /**
   * @brief 围绕名义状态和控制轨迹构造离散 LQ 近似。
   *
   * 中间节点由 Transcription 层根据连续/离散动力学模式生成；终端节点仅包含
   * 终端目标函数和终端约束的二次近似。
   *
   * 该方法会更新以下变量：
   * - 线性化后的离散动力学模型
   * - 二次化目标函数
   * - 增广拉格朗日项产生的约束惩罚近似
   */
  void approximateOptimalControlProblem() {
    approximateIntermediateLQ(nominalDualData_.dualSolution,
                              nominalPrimalData_);

    /*
     * 计算终端时刻的 Heuristics 函数，并调用 shiftHessian
     * 处理 Heuristics 的二阶导数。
     */
    ModelData_t& modelData = nominalPrimalData_.modelDataFinalTime;
    const Scalar& time =
        nominalPrimalData_.primalSolution.timeTrajectory_.back();
    const StateVector_t& state =
        nominalPrimalData_.primalSolution.stateTrajectory_.back();
    const FinalMultiplierCollection_t& multiplier =
        nominalDualData_.dualSolution.final;
    modelData = LinearQuadraticApproximator_t::approximateFinalLQ(
        optimalControlProblem_, targetTrajectory_, time, state, multiplier);

    // 修正终端时刻 Hessian。
    shiftHessian(ddpSettings_.lineSearch.hessianCorrectionStrategy,
                 modelData.cost.dfdxx,
                 ddpSettings_.lineSearch.hessianCorrectionMultiple);
  }

  /**
   * @brief 计算所有中间节点处最优控制问题的离散 LQ 近似。
   *
   * @param [in] dualSolution 对偶解。
   * @param [in,out] primalData 原始数据（读轨迹，写 modelDataTrajectory）。
   */
  void approximateIntermediateLQ(const DualSolution_t& dualSolution,
                                 PrimalDataContainer_t& primalData) {
    // 创建别名。
    const TimeTrajectory_t& timeTrajectory =
        primalData.primalSolution.timeTrajectory_;
    const StateTrajectory_t& stateTrajectory =
        primalData.primalSolution.stateTrajectory_;
    const InputTrajectory_t& inputTrajectory =
        primalData.primalSolution.inputTrajectory_;
    const IntermediateMultiplierTrajectory_t& multiplierTrajectory =
        dualSolution.intermediates;
    ModelDataTrajectory_t& modelDataTrajectory = primalData.modelDataTrajectory;

    for (size_t timeIndex = 0; timeIndex < PredictLength; ++timeIndex) {
      const Scalar timeStep =
          timeTrajectory[timeIndex + 1] - timeTrajectory[timeIndex];
      transcription_.approximateIntermediateLQ(
          optimalControlProblem_, targetTrajectory_, timeTrajectory[timeIndex],
          stateTrajectory[timeIndex], inputTrajectory[timeIndex], timeStep,
          multiplierTrajectory[timeIndex], modelDataTrajectory[timeIndex]);
    };
  }

  /**
   * @brief 计算哈密顿量对控制的 Hessian（Hm = dfduu + B^T Sm
   * B），并可能被搜索策略修正。
   * @param [in] modelData 当前节点模型数据。
   * @param [in] Sm 当前价值函数二阶导（Riccati 矩阵）。
   * @return 哈密顿量 Hessian 矩阵 Hm。
   */
  HmMatrix_t computeHamiltonianHessian(const ModelData_t& modelData,
                                       const SmMatrix_t& Sm) const {
    const Matrix<Scalar, UDim, XDim> BmTransSm =
        modelData.dynamics.dfdu.transpose() * Sm;
    HmMatrix_t Hm = modelData.cost.dfduu;
    Hm += BmTransSm * modelData.dynamics.dfdu;
    return lineSearchStrategy_.augmentHamiltonianHessian(modelData, Hm);
  }

  /**
   * @brief 计算 Riccati 递推所需的 Hamiltonian Hessian、LLT 分解和
   * Riccati 方程修正项。
   *
   * @param [in] modelData 当前节点模型数据。
   * @param [in] Sm 下一时刻 Riccati 矩阵。
   * @param [out] riccatiModification Riccati 修正项。
   */
  void prepareRiccatiModification(
      const ModelData_t& modelData, const SmMatrix_t& Sm,
      RiccatiModification_t& riccatiModification) const {
    // 计算 Hamiltonian 的 Hessian。
    riccatiModification.time_ = modelData.time;
    riccatiModification.hamiltonianHessian_ =
        computeHamiltonianHessian(modelData, Sm);

    // 计算Hessian的LLT分解
    riccatiModification.HmLLT_.Decomposition(
        riccatiModification.hamiltonianHessian_);

    // 计算 deltaQm、deltaGv、deltaGm。
    lineSearchStrategy_.computeRiccatiModification(
        modelData, riccatiModification.deltaQm_);
  }

  /**
   * @brief 从终端到初始时刻顺序求解 Riccati 方程，得到 value function 与
   * Lv/Km。
   * @param [in] finalValueFunction 终端价值函数二次近似（Sm=dfdxx, Sv=dfdx,
   * s=f）。
   * @param [out] controllerGainTrajectory 控制器反馈增益数组；仅写入
   * 0..PredictLength-1，终端点由 calculateController() 补齐。
   * @param [out] controllerDeltaBiasTrajectory 控制器前馈更新数组；仅写入
   * 0..PredictLength-1，终端点由 calculateController() 补齐。
   */
  void solveSequentialRiccatiEquations(
      const ValueFunctionQuadraticApproximation_t& finalValueFunction,
      ControllerGainTrajectory_t& controllerGainTrajectory,
      ControllerDeltaBiasTrajectory_t& controllerDeltaBiasTrajectory) {
    nominalDualData_.valueFunctionTrajectory.back() = finalValueFunction;

    riccatiEquationsWorker(finalValueFunction, controllerGainTrajectory,
                           controllerDeltaBiasTrajectory);
  }

  /**
   * 求解给定索引分区中的 Riccati 方程和 type_1 约束误差校正
   * 补偿。
   *
   * @param [in] finalValueFunction 终端价值函数，用于反向递推的初值。
   * @param [out] controllerGainTrajectory 控制器反馈增益数组。
   * @param [out] controllerDeltaBiasTrajectory 控制器前馈更新数组。
   */
  void riccatiEquationsWorker(
      const ValueFunctionQuadraticApproximation_t& finalValueFunction,
      ControllerGainTrajectory_t& controllerGainTrajectory,
      ControllerDeltaBiasTrajectory_t& controllerDeltaBiasTrajectory) {
    /*
     * 求解 Riccati 方程。
     */
    const ValueFunctionQuadraticApproximation_t* valueFunctionNext =
        &finalValueFunction;

    // Riccati 只产生区间控制律 K_0..K_{N-1} 和 d_0..d_{N-1}；
    // controller 的末端点用于时间插值，稍后由倒数第二个控制律补齐。
    int curIndex = PredictLength - 1;
    constexpr int stopIndex = 0;
    while (curIndex >= stopIndex) {
      LvVector_t& curLv = controllerDeltaBiasTrajectory[curIndex];
      KmMatrix_t& curKm = controllerGainTrajectory[curIndex];
      RiccatiModification_t& curRiccatiModification =
          nominalDualData_.riccatiModificationTrajectory[curIndex];
      const ModelData_t& curModelData =
          nominalPrimalData_.modelDataTrajectory[curIndex];

      SmMatrix_t& curSm =
          nominalDualData_.valueFunctionTrajectory[curIndex].dfdxx;
      SvVector_t& curSv =
          nominalDualData_.valueFunctionTrajectory[curIndex].dfdx;
      Scalar& curs = nominalDualData_.valueFunctionTrajectory[curIndex].f;

      prepareRiccatiModification(curModelData, valueFunctionNext->dfdxx,
                                 curRiccatiModification);

      riccatiEquationsSolver_.computeMap(
          curModelData, curRiccatiModification, valueFunctionNext->dfdxx,
          valueFunctionNext->dfdx, valueFunctionNext->f, curKm, curLv, curSm,
          curSv, curs);
      valueFunctionNext = &(nominalDualData_.valueFunctionTrajectory[curIndex]);

      --curIndex;
    }  // while 循环。
  }

  /**
   * 计算控制器。该方法使用以下变量，并会
   * 修改 unoptimizedController_。
   */
  void calculateController() {
    unoptimizedController_.timeStamp_ =
        nominalPrimalData_.primalSolution.timeTrajectory_;
    optimizedPrimalSolution_.controller_.timeStamp_ =
        unoptimizedController_.timeStamp_;

    for (size_t timeIndex = 0; timeIndex < PredictLength; ++timeIndex) {
      calculateControllerWorker(timeIndex, nominalPrimalData_,
                                unoptimizedController_);
    }

    // 由于最后一个时间戳的控制器无效，如果最后时刻
    // 不是事件时刻，则使用倒数第二个时刻的控制策略
    // 作为最后时刻控制策略；finalTimeIsNotAnEvent && 至少有两个
    // 时间戳。
    if (unoptimizedController_.size() >= 2) {
      const size_t secondToLastIndex = unoptimizedController_.size() - 2u;
      unoptimizedController_.gainArray_.back() =
          unoptimizedController_.gainArray_[secondToLastIndex];
      unoptimizedController_.biasArray_.back() =
          unoptimizedController_.biasArray_[secondToLastIndex];
      unoptimizedController_.deltaBiasArray_.back() =
          unoptimizedController_.deltaBiasArray_[secondToLastIndex];
    }
  }

  /**
   * 使用原始解计算 timeIndex 对应的控制器，并写回 dstController。
   *
   * @param [in] timeIndex 当前时间索引。
   * @param [in] primalData 用于计算控制器的原始数据。
   * @param [out] dstController 输出的控制器（增益、偏置、deltaBias
   * 写入对应索引）。
   */
  void calculateControllerWorker(size_t timeIndex,
                                 const PrimalDataContainer_t& primalData,
                                 LinearController_t& dstController) {
    const StateVector_t& nominalState =
        primalData.primalSolution.stateTrajectory_[timeIndex];
    const InputVector_t& nominalInput =
        primalData.primalSolution.inputTrajectory_[timeIndex];

    // gainArray_/deltaBiasArray_ 已由 Riccati 递推写入；这里仅把
    // u = Kx + uff 转换成经过名义点的 affine bias：uff = u_nominal -
    // Kx_nominal。
    dstController.biasArray_[timeIndex] = nominalInput;
    dstController.biasArray_[timeIndex] -=
        dstController.gainArray_[timeIndex] * nominalState;
  }

  /** 基于当前 LQ 解更新优化后的原始解和对偶
   * 解。 */
  void takePrimalDualStep(Scalar lqModelExpectedCost) {
    // 更新原始解：运行搜索策略并找到最优 stepLength。
    SearchStrategySolutionRef_t solution(
        optimizedDualSolution_, optimizedPrimalSolution_,
        optimizedProblemMetrics_, performanceIndex_);
    const bool success = lineSearchStrategy_.run(
        {initTime_, finalTime_}, initState_, lqModelExpectedCost,
        unoptimizedController_, nominalDualData_.dualSolution, solution);

    // 更新对偶解。
    if (success) {
      DualSolutionRef_t DualSolutionRef = optimizedDualSolution_;
      updateDualSolution(optimalControlProblem_, optimizedPrimalSolution_,
                         optimizedProblemMetrics_, DualSolutionRef);
      performanceIndex_ = computeRolloutPerformanceIndex(
          optimizedPrimalSolution_.timeTrajectory_, optimizedProblemMetrics_);
      performanceIndex_.merit = calculateRolloutMerit(performanceIndex_);
    }

    // 如果失败，则使用名义解；为保持缓存数据一致性，所有
    // 缓存都应保持不变。
    if (!success) {
      optimizedDualSolution_ = nominalDualData_.dualSolution;
      optimizedPrimalSolution_ = nominalPrimalData_.primalSolution;
      optimizedProblemMetrics_ = nominalPrimalData_.problemMetrics;
      performanceIndex_ = performanceIndexLast_;
    }
  }

  ValueFunctionQuadraticApproximation_t getValueFunctionImpl(
      const Scalar time, const StateVector_t& state,
      const PrimalSolution_t& primalSolution,
      const ValueFunctionTrajectory_t& valueFunctionTrajectory) const {
    // 结果。
    ValueFunctionQuadraticApproximation_t valueFunction;
    valueFunction.dfdu.setZero();
    valueFunction.dfdux.setZero();
    valueFunction.dfduu.setZero();
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, primalSolution.timeTrajectory_);
    valueFunction.f = LinearInterpolation::interpolate(
        indexAlpha, valueFunctionTrajectory,
        +[](const ValueFunctionTrajectory_t& trajectory,
            size_t ind) -> const Scalar& { return trajectory[ind].f; });
    valueFunction.dfdx = LinearInterpolation::interpolate(
        indexAlpha, valueFunctionTrajectory,
        +[](const ValueFunctionTrajectory_t& trajectory, size_t ind)
            -> const StateVector_t& { return trajectory[ind].dfdx; });
    valueFunction.dfdxx = LinearInterpolation::interpolate(
        indexAlpha, valueFunctionTrajectory,
        +[](const ValueFunctionTrajectory_t& trajectory,
            size_t ind) -> const SmMatrix_t& { return trajectory[ind].dfdxx; });

    // 围绕查询状态重新居中。
    const StateVector_t xNominal = LinearInterpolation::interpolate(
        indexAlpha, primalSolution.stateTrajectory_);
    const StateVector_t deltaX = state - xNominal;
    const StateVector_t SmDeltaX = valueFunction.dfdxx * deltaX;
    valueFunction.f += deltaX.dot(Scalar(0.5) * SmDeltaX + valueFunction.dfdx);
    valueFunction.dfdx += SmDeltaX;  // 在更新 f 后调整 dfdx！

    return valueFunction;
  }

  ValueFunctionQuadraticApproximation_t getValueFunctionImpl(
      size_t timeIndex, const StateVector_t& state,
      const PrimalSolution_t& primalSolution,
      const ValueFunctionTrajectory_t& valueFunctionTrajectory) const {
    assert(timeIndex <= PredictLength);

    ValueFunctionQuadraticApproximation_t valueFunction =
        valueFunctionTrajectory[timeIndex];
    valueFunction.dfdu.setZero();
    valueFunction.dfdux.setZero();
    valueFunction.dfduu.setZero();

    const StateVector_t deltaX =
        state - primalSolution.stateTrajectory_[timeIndex];
    const StateVector_t SmDeltaX = valueFunction.dfdxx * deltaX;
    valueFunction.f += deltaX.dot(Scalar(0.5) * SmDeltaX + valueFunction.dfdx);
    valueFunction.dfdx += SmDeltaX;  // 在更新 f 后调整 dfdx！

    return valueFunction;
  }

  DDPSettings<Scalar> ddpSettings_{};

  // 最优控制问题定义由外部持有，本求解器仅保存引用。
  OptimalControlProblem_t& optimalControlProblem_;

  Rollout_t rollout_;

  InitializerRollout_t initializerRollout_;

  LineSearchStrategy_t lineSearchStrategy_;

  Transcription_t transcription_;

  TargetTrajectories_t targetTrajectory_;

  // 当前求解时间区间和初始条件。
  Scalar initTime_{Scalar(0.0)};
  Scalar finalTime_{Scalar(0.0)};
  Scalar lastFinalTime_{Scalar(0.0)};
  StateVector_t initState_;

  // 当前反向递推使用的名义数据。
  PrimalDataContainer_t nominalPrimalData_;
  DualDataContainer_t nominalDualData_;

  // 当前接受的最佳解及其 rollout 指标。
  PrimalSolution_t optimizedPrimalSolution_;
  DualSolution_t optimizedDualSolution_;
  ProblemMetrics_t optimizedProblemMetrics_;

  // 线搜索前由 LQ 解计算出的控制器。
  LinearController_t unoptimizedController_;

  // Riccati 递推工作区。
  DiscreteTimeRiccatiEquations_t riccatiEquationsSolver_;

  // 性能和迭代记录。
  PerformanceIndex_t performanceIndex_;
  PerformanceIndex_t performanceIndexLast_;
  size_t totalNumIterations_{0};
  size_t totalNumRuns_{0};
};