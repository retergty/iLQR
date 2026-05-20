/**
 * @file iLQR.hpp
 * @brief 离散时间 iLQR 求解器：在名义轨迹上做 LQ 近似并迭代求解
 * Riccati，支持线搜索与约束投影。
 */
#pragma once
#include "ChangeOfInputVariables.hpp"
#include "ControlledSystemBase.hpp"
#include "Controller.hpp"
#include "Cost.hpp"
#include "DDPData.hpp"
#include "DDPSetting.hpp"
#include "DefaultInitializer.hpp"
#include "DiscreteTimeRiccatiEquations.hpp"
#include "DualSolution.hpp"
#include "HessianCorrection.hpp"
#include "Initializer.hpp"
#include "InitializerRollout.hpp"
#include "LineSearchStrategy.hpp"
#include "LinearAlgebra.hpp"
#include "LinearController.hpp"
#include "LinearQuadraticApproximator.hpp"
#include "OptimalControlProblem.hpp"
#include "OptimalControlProblemHelperFunction.hpp"
#include "PerformanceIndex.hpp"
#include "PrimalSolution.hpp"
#include "ProblemMetrics.hpp"
#include "RolloutBase.hpp"
#include "SearchStrategyBase.hpp"
#include "SensitivityIntegrator.hpp"
#include "TimeTriggeredRollout.hpp"
#include "TrapezoidalIntegration.hpp"
#include "Types.hpp"
#include "iLQRDescriptor.hpp"
#include "iLQRTypes.hpp"

/**
 * @brief 迭代线性二次调节器（iLQR）：基于名义轨迹的 LQ 近似与离散时间 Riccati
 * 反向递推的 DDP 求解器。
 * @tparam Scalar 标量类型（如 double）。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 * @tparam PredictLength 预测步数（时间节点数为 PredictLength+1）。
 * @tparam StateEqConstrains 状态等式约束数（中间时刻）。
 * @tparam StateIneqConstrains 状态不等式约束数（中间时刻）。
 * @tparam StateInputEqConstrains 状态-输入等式约束数。
 * @tparam StateInputIneqConstrains 状态-输入不等式约束数。
 * @tparam FinalStateEqConstrains 终端状态等式约束数。
 * @tparam FinalStateIneqConstrains 终端状态不等式约束数。
 */
template <typename Descriptor>
class iLQR {
 public:
  using Types = iLQRTypes<Descriptor>;
  using Scalar = typename Types::Scalar;

  static constexpr int XDim = Types::XDim;
  static constexpr int UDim = Types::UDim;
  static constexpr std::size_t PredictLength = Types::PredictLength;
  static constexpr int StateEqConstrains = Types::StateEq;
  static constexpr int StateIneqConstrains = Types::StateIneq;
  static constexpr int StateInputEqConstrains = Types::StateInputEq;
  static constexpr int StateInputIneqConstrains = Types::StateInputIneq;
  static constexpr int FinalStateEqConstrains = Types::FinalStateEq;
  static constexpr int FinalStateIneqConstrains = Types::FinalStateIneq;

  using OptimalControlProblem_t = typename Types::OptimalControlProblem_t;
  using ControlledSystemBase_t = typename Types::ControlledSystemBase_t;
  using SystemDynamicsBase_t = typename Types::SystemDynamicsBase_t;

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

  using TimeTriggeredRollout_t = typename Types::TimeTriggeredRollout_t;
  using RolloutTrajectoryPointer_t = typename Types::RolloutTrajectoryPointer_t;

  using SearchStrategySolution_t = typename Types::SearchStrategySolution_t;
  using SearchStrategySolutionRef_t =
      typename Types::SearchStrategySolutionRef_t;

  using LinearController_t = typename Types::LinearController_t;
  using SearchStrategyBase_t = typename Types::SearchStrategyBase_t;
  using LineSearchStrategy_t = LineSearchStrategy<Descriptor>;
  using RiccatiModification_t = typename Types::RiccatiModification_t;

  using TimeTrajectory_t = typename Types::TimeTrajectory_t;
  using StateTrajectory_t = typename Types::StateTrajectory_t;
  using InputTrajectory_t = typename Types::InputTrajectory_t;
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
  using IntermediatePerformanceIndexTrajectory_t =
      typename Types::IntermediatePerformanceIndexTrajectory_t;

  using EK2DynamicsDiscretizer_t = typename Types::EK2DynamicsDiscretizer_t;

  using ValueFunctionQuadraticApproximation_t =
      typename Types::ValueFunctionQuadraticApproximation_t;
  using ValueFunctionTrajectory_t = typename Types::ValueFunctionTrajectory_t;
  using DiscreteTimeRiccatiEquations_t =
      typename Types::DiscreteTimeRiccatiEquations_t;

  /**
   * @brief 构造 iLQR 求解器，绑定动力学与初始化器。
   * @param [in] systemPtr 系统动力学指针（用于 rollout 与离散化），不可为
   * nullptr。
   * @param [in] initializer 轨迹初始化器，用于填补无控制器的时间段。
   */
  iLQR(const DDPSettings<Scalar>& ddp_setting, SystemDynamicsBase_t* systemPtr,
       Initializer_t* initializer)
      : ddpSettings_(ddp_setting),
        rollout_(systemPtr, ddpSettings_.timeStep_),
        initializerRollout_(*initializer, ddpSettings_.timeStep_),
        lineSearchStrategy_(makeLineSearchSettings(ddp_setting), *this) {
    optimalControlProblem_.dynamicsPtr = systemPtr;
    // set zero solution
    optimizedPrimalSolution_.clear();
    optimizedDualSolution_.clear();
  };

  /**
   * The main routine of solver which runs the optimizer for a given initial
   * state, initial time, and final time.
   *
   * @param [in] initTime 初始时间。
   * @param [in] initState 初始状态。
   */
  void run(Scalar initTime, const StateVector_t& initState) {
    // initialize parameters
    initTime_ = initTime;
    initState_ = initState;
    lastFinalTime_ = finalTime_;
    finalTime_ = initTime + ddpSettings_.timeStep_ * (PredictLength);
    const size_t initIteration = totalNumIterations_;

    // optimized --> nominal: initializes the nominal primal and dual solutions
    // based on the optimized ones
    bool initialSolutionExists =
        initializePrimalSolution();  // true if the rollout is not purely from
                                     // the Initializer
    initializeDualSolutionAndMetrics();

    performanceIndexLast_ = performanceIndex_;
    // convergence variables of the main loop
    bool isConverged = false;

    // DDP main loop
    while (true) {
      // nominal --> nominal: constructs the LQ problem around the nominal
      // trajectories
      approximateOptimalControlProblem();

      // nominal --> nominal: solves the LQ problem
      avgTimeStepBP_ = solveSequentialRiccatiEquations(
          nominalPrimalData_.modelDataFinalTime.cost);

      // calculate controller and store the result in unoptimizedController_
      calculateController();

      // the expected cost/merit calculated by the Riccati solution is not
      // reliable
      const auto lqModelExpectedCost =
          initialSolutionExists
              ? nominalDualData_.valueFunctionTrajectory.front().f
              : performanceIndex_.merit;

      // nominal --> optimized: based on the current LQ solution updates the
      // optimized primal and dual solutions
      takePrimalDualStep(lqModelExpectedCost);

      // iteration info
      ++totalNumIterations_;

      // check convergence
      isConverged = lineSearchStrategy_.checkConvergence(
          !initialSolutionExists, performanceIndexLast_, performanceIndex_);
      performanceIndexLast_ = performanceIndex_;
      initialSolutionExists = true;

      if (isConverged || (totalNumIterations_ - initIteration) ==
                             ddpSettings_.maxNumIterations_) {
        break;
      } else {
        // optimized --> nominal: use the optimized solution as the nominal for
        // the next iteration
        optimizedDualSolution_.swap(nominalDualData_.dualSolution);
        optimizedPrimalSolution_.swap(nominalPrimalData_.primalSolution);
        optimizedProblemMetrics_.swap(nominalPrimalData_.problemMetrics);
      }
    }  // end of while loop
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
    optimalControlProblem_.timeTrajectory = timeTrajectory;
    optimalControlProblem_.stateTrajectory = stateTrajectory;
    optimalControlProblem_.inputTrajectory = inputTrajectory;
  }

  const PrimalSolution_t& primalSolution() const {
    return optimizedPrimalSolution_;
  }
  const PerformanceIndex_t performanceIndex() const {
    return performanceIndex_;
  }

 private:
  static LineSearchSettings<Scalar> makeLineSearchSettings(
      const DDPSettings<Scalar>& ddpSettings) {
    LineSearchSettings<Scalar> settings = ddpSettings.lineSearch_;
    settings.minRelCost = ddpSettings.minRelCost_;
    return settings;
  }

  /** Initializes the nominal primal based on the optimized ones.
   * @return True if the rollout is not purely from the Initializer.
   */
  bool initializePrimalSolution() {
    // try to initialize with controller
    int numSteps = 0;
    bool ret = false;
    if (lastFinalTime_ > initTime_) {
      numSteps = rolloutInitialController(optimizedPrimalSolution_,
                                          nominalPrimalData_.primalSolution);
      ret = true;
    }
    // past too long before last run.
    rolloutInitializer(nominalPrimalData_.primalSolution, numSteps);
    return ret;
  }

  /**
   * Forward integrate the system dynamics with the controller in
   * inputPrimalSolution. In general, it uses the given control policies and the
   * initial state, to integrate the system dynamics in the time period
   * [initTime, finalTime]. However, if inputPrimalSolution's controller does
   * not cover the period [initTime, finalTime], it will use the controller till
   * the final time of the controller
   *
   * @param [in] inputPrimalSolution 其控制器将用于前向积分。
   * @param [out] outputPrimalSolution 得到的原始解（轨迹与控制器）。
   * @return 被覆盖的步数。
   */
  int rolloutInitialController(PrimalSolution_t& inputPrimalSolution,
                               PrimalSolution_t& outputPrimalSolution) {
    // Ensure that finalTime is included by adding a fraction of dt such that: N
    // * dt <= finalTime < (N + 1) * dt.
    Scalar finalTimeLocal = std::min(lastFinalTime_, finalTime_) +
                            static_cast<Scalar>(0.01) * ddpSettings_.timeStep_;
    int numSteps = std::min(
        static_cast<int>((finalTimeLocal - initTime_) / ddpSettings_.timeStep_),
        static_cast<int>(PredictLength));

    outputPrimalSolution.controller_ = inputPrimalSolution.controller_;
    rolloutTrajectory(rollout_, initTime_, initState_, lastFinalTime_,
                      outputPrimalSolution);

    return numSteps;
  }

  /**
   * It will check the content of the primalSolution, and if its final time is
   * smaller than the current solver finalTime_, it will concatenate it with the
   * result of Initializer.
   */
  void rolloutInitializer(PrimalSolution_t& primalSolution, int numSteps) {
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
   * Initializes the nominal dual solutions based on the optimized ones and
   * nominal primal solution. Moreover, it updates ProblemMetrics.
   */
  void initializeDualSolutionAndMetrics() {
    // initialize dual solution
    initializeDualSolution(
        optimalControlProblem_, nominalPrimalData_.primalSolution,
        optimizedDualSolution_, nominalDualData_.dualSolution);

    computeRolloutMetrics(
        optimalControlProblem_, nominalPrimalData_.primalSolution,
        nominalDualData_.dualSolution, nominalPrimalData_.problemMetrics);

    // calculates rollout merit
    performanceIndex_ = computeRolloutPerformanceIndex(
        nominalPrimalData_.primalSolution.timeTrajectory_,
        nominalPrimalData_.problemMetrics);
    performanceIndex_.merit = calculateRolloutMerit(performanceIndex_);
  }

  /**
   * Approximates the nonlinear problem as a linear-quadratic problem around the
   * nominal state and control trajectories. This method updates the following
   * variables:
   * 	- linearized system model and constraints
   * 	- quadratized cost function
   * 	- as well as the constrained coefficients of
   * 		- linearized system model
   * 		- quadratized intermediate cost function
   * 		- quadratized final cost
   */
  void approximateOptimalControlProblem() {
    approximateIntermediateLQ(nominalDualData_.dualSolution,
                              nominalPrimalData_);

    /*
     * compute the Heuristics function at the final time. Also call shiftHessian
     * on the Heuristics 2nd order derivative.
     */
    ModelData_t& modelData = nominalPrimalData_.modelDataFinalTime;
    const Scalar& time =
        nominalPrimalData_.primalSolution.timeTrajectory_.back();
    const StateVector_t& state =
        nominalPrimalData_.primalSolution.stateTrajectory_.back();
    const FinalMultiplierCollection_t& multiplier =
        nominalDualData_.dualSolution.final;
    modelData = approximator_.approximateFinalLQ(optimalControlProblem_, time,
                                                 state, multiplier);

    // shift Hessian for final time
    if (ddpSettings_.strategy_ == SearchStrategyType::LINE_SEARCH) {
      shiftHessian(ddpSettings_.lineSearch_.hessianCorrectionStrategy,
                   modelData.cost.dfdxx,
                   ddpSettings_.lineSearch_.hessianCorrectionMultiple);
    }
  }

  /**
   * Calculates an LQ approximate of the optimal control problem for the nodes.
   *
   * @param [in] dualSolution 对偶解。
   * @param [in,out] primalData 原始数据（读轨迹，写 modelDataTrajectory）。
   */
  void approximateIntermediateLQ(const DualSolution_t& dualSolution,
                                 PrimalDataContainer_t& primalData) {
    // create alias
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
      // approximate continuous LQ for the given time index
      ModelData_t continuousTimeModelData =
          approximator_.approximateIntermediateLQ(
              optimalControlProblem_, timeTrajectory[timeIndex],
              stateTrajectory[timeIndex], inputTrajectory[timeIndex],
              multiplierTrajectory[timeIndex]);

      // TO DO:checking the numerical properties

      // // discretize LQ problem
      // const Scalar timeStep = (timeIndex + 1 < timeTrajectory.size()) ?
      // (timeTrajectory[timeIndex + 1] - timeTrajectory[timeIndex]) : 0.0; if
      // (!numerics::almost_eq(timeStep, 0.0))
      // {
      //     discreteLQWorker(*optimalControlProblem_.dynamicsPtr,
      //     timeTrajectory[timeIndex], stateTrajectory[timeIndex],
      //                      inputTrajectory[timeIndex], timeStep,
      //                      continuousTimeModelData,
      //                      modelDataTrajectory[timeIndex]);
      // }
      // else
      // {
      //     modelDataTrajectory[timeIndex] = continuousTimeModelData;
      // }

      const Scalar timeStep =
          timeTrajectory[timeIndex + 1] - timeTrajectory[timeIndex];
      discreteLQWorker(*optimalControlProblem_.dynamicsPtr,
                       timeTrajectory[timeIndex], stateTrajectory[timeIndex],
                       inputTrajectory[timeIndex], timeStep,
                       continuousTimeModelData, modelDataTrajectory[timeIndex]);
    };
  }

  /**
   * Calculates the discrete-time LQ approximation from the continuous-time LQ
   * approximation.
   *
   * @param [in] system 系统动力学。
   * @param [in] time 当前时刻 t_k。
   * @param [in] state 当前状态 x_k。
   * @param [in] input 当前输入 u_k。
   * @param [in] timeStep x_k 到 x_{k+1} 的时间步长。
   * @param [in] continuousTimeModelData 连续时间 LQ 模型数据。
   * @param [out] modelData 离散化后的模型数据。
   */
  void discreteLQWorker(SystemDynamicsBase<Scalar, XDim, UDim>& system,
                        Scalar time, const StateVector_t& state,
                        const InputVector_t& input, Scalar timeStep,
                        const ModelData_t& continuousTimeModelData,
                        ModelData_t& modelData) {
    modelData.time = continuousTimeModelData.time;

    // linearize system dynamics
    modelData.dynamics = discretizer_.sensitivityDiscretize(system, time, state,
                                                            input, timeStep);
    modelData.dynamics.f.setZero();  // why?

    // quadratic approximation to the cost function
    modelData.cost = continuousTimeModelData.cost * timeStep;
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
   * @brief 由 Hm 的 UUT 分解得到 inv(Hm)，用作约束零空间投影（无约束时即
   * inv(Hm)）。
   * @param [in] Hm 哈密顿量对控制的 Hessian（正定）。
   * @param [out] constraintNullProjector 输出投影矩阵（当前实现为 inv(Hm) 的
   * UUT 因子）。
   */
  void computeProjections(const HmMatrix_t& Hm,
                          HmMatrix_t& constraintNullProjector) const {
    // UUT decomposition of inv(Hm)
    HmMatrix_t HmInvUmUmT;

    LinearAlgebra::computeInverseMatrixUUT(Hm, HmInvUmUmT);

    // compute DmDagger, DmDaggerTHmDmDaggerUUT, HmInverseConstrainedLowRank
    constraintNullProjector = HmInvUmUmT;
  }

  /**
   * @brief 将 LQ 模型按 u = Pu * tilde{u}
   * 做变量替换，得到投影后的动力学与代价系数。
   * @param [in] modelData 原始模型数据。
   * @param [in] constraintNullProjector 投影矩阵 Pu（约束零空间）。
   * @param [out] projectedModelData 投影后的模型数据。
   */
  void projectLQ(const ModelData_t& modelData,
                 const HmMatrix_t& constraintNullProjector,
                 ModelData_t& projectedModelData) const {
    // dimensions and time
    projectedModelData.time = modelData.time;

    // Change of variables u = Pu * tilde{u}
    // Pu = constraintNullProjector;

    // dynamics
    projectedModelData.dynamics = modelData.dynamics;
    changeOfInputVariables(projectedModelData.dynamics,
                           constraintNullProjector);

    // cost
    projectedModelData.cost = modelData.cost;
    changeOfInputVariables(projectedModelData.cost, constraintNullProjector);
  }

  /**
   * Takes the following steps: (1) Computes the Hessian of the Hamiltonian
   * (i.e., Hm) (2) Based on Hm, it calculates the range space and the null
   * space projections of the input-state equality constraints. (3) Based on
   * these two projections, defines the projected LQ model. (4) Finally, defines
   * the Riccati equation modifiers based on the search strategy.
   *
   * @param [in] modelData 当前节点模型数据。
   * @param [in] Sm 下一时刻 Riccati 矩阵。
   * @param [out] projectedModelData 投影后的模型数据。
   * @param [out] riccatiModification Riccati 修正项。
   */
  void computeProjectionAndRiccatiModification(
      const ModelData_t& modelData, const SmMatrix_t& Sm,
      ModelData_t& projectedModelData,
      RiccatiModification_t& riccatiModification) const {
    // compute the Hamiltonian's Hessian
    riccatiModification.time_ = modelData.time;
    riccatiModification.hamiltonianHessian_ =
        computeHamiltonianHessian(modelData, Sm);

    // compute projectors
    computeProjections(riccatiModification.hamiltonianHessian_,
                       riccatiModification.constraintNullProjector_);

    // project LQ
    projectLQ(modelData, riccatiModification.constraintNullProjector_,
              projectedModelData);

    // compute deltaQm, deltaGv, deltaGm
    lineSearchStrategy_.computeRiccatiModification(
        projectedModelData, riccatiModification.deltaQm_);
  }

  /**
   * @brief 从终端到初始时刻顺序求解 Riccati 方程，得到 value function 与
   * projected Lv/Km。
   * @param [in] finalValueFunction 终端价值函数二次近似（Sm=dfdxx, Sv=dfdx,
   * s=f）。
   * @return 平均时间步长（用于统计）。
   */
  Scalar solveSequentialRiccatiEquations(
      const ValueFunctionQuadraticApproximation_t& finalValueFunction) {
    const ModelData_t& finalModelData = nominalPrimalData_.modelDataFinalTime;
    RiccatiModification_t& finalRiccatiModification =
        nominalDualData_.riccatiModificationTrajectory.back();
    ModelData_t& finalProjectedModelData =
        nominalDualData_.projectedModelDataTrajectory.back();
    LvVector_t& finalProjectedLvFinal = projectedLvTrajectoryStock_.back();
    KmMatrix_t& finalProjectedKmFinal = projectedKmTrajectoryStock_.back();

    SmMatrix_t SmDummy;
    SmDummy.setZero();

    computeProjectionAndRiccatiModification(finalModelData, SmDummy,
                                            finalProjectedModelData,
                                            finalRiccatiModification);

    // projected feedforward
    finalProjectedLvFinal = -finalProjectedModelData.cost.dfdu;
    // last
    // finalProjectedLvFinal -=
    // finalProjectedModelData.dynamics.dfdu.transpose() *
    // finalValueFunction.dfdx;

    // projected feedback
    finalProjectedKmFinal = -finalProjectedModelData.cost.dfdux;
    // finalProjectedKmFinal -=
    // finalProjectedModelData.dynamics.dfdu.transpose() *
    // finalValueFunction.dfdxx;

    return solveSequentialRiccatiEquationsImpl(finalValueFunction);
  }

  /**
   * The implementation for solving Riccati equations for all the partitions.
   *
   * @param [in] finalValueFunction 终端价值函数（Sm=dfdxx, Sv=dfdx, s=f）。
   * @return 平均时间步长。
   */
  Scalar solveSequentialRiccatiEquationsImpl(
      const ValueFunctionQuadraticApproximation_t& finalValueFunction) {
    nominalDualData_.valueFunctionTrajectory.back() = finalValueFunction;

    riccatiEquationsWorker(finalValueFunction);

    // average time step
    return (finalTime_ - initTime_) / PredictLength;
  }

  /**
   * Solves a Riccati equations and type_1 constraints error correction
   * compensation for the partition in the given index.
   *
   * @param [in] finalValueFunction 终端价值函数，用于反向递推的初值。
   */
  void riccatiEquationsWorker(
      const ValueFunctionQuadraticApproximation_t& finalValueFunction) {
    /*
     * solving the Riccati equations
     */
    const ValueFunctionQuadraticApproximation_t* valueFunctionNext =
        &finalValueFunction;

    int curIndex = PredictLength - 1;
    constexpr int stopIndex = 0;
    while (curIndex >= stopIndex) {
      LvVector_t& curProjectedLv = projectedLvTrajectoryStock_[curIndex];
      KmMatrix_t& curProjectedKm = projectedKmTrajectoryStock_[curIndex];
      ModelData_t& curProjectedModelData =
          nominalDualData_.projectedModelDataTrajectory[curIndex];
      RiccatiModification_t& curRiccatiModification =
          nominalDualData_.riccatiModificationTrajectory[curIndex];
      const ModelData_t& curModelData =
          nominalPrimalData_.modelDataTrajectory[curIndex];

      SmMatrix_t& curSm =
          nominalDualData_.valueFunctionTrajectory[curIndex].dfdxx;
      SvVector_t& curSv =
          nominalDualData_.valueFunctionTrajectory[curIndex].dfdx;
      Scalar& curs = nominalDualData_.valueFunctionTrajectory[curIndex].f;

      computeProjectionAndRiccatiModification(
          curModelData, valueFunctionNext->dfdxx, curProjectedModelData,
          curRiccatiModification);

      riccatiEquationsSolver_.computeMap(
          curProjectedModelData, curRiccatiModification,
          valueFunctionNext->dfdxx, valueFunctionNext->dfdx,
          valueFunctionNext->f, curProjectedKm, curProjectedLv, curSm, curSv,
          curs);
      valueFunctionNext = &(nominalDualData_.valueFunctionTrajectory[curIndex]);

      --curIndex;
    }  // while
  }

  /**
   * Calculates the controller. This method uses the following variables. The
   * method modifies unoptimizedController_.
   */
  void calculateController() {
    unoptimizedController_.timeStamp_ =
        nominalPrimalData_.primalSolution.timeTrajectory_;
    optimizedPrimalSolution_.controller_.timeStamp_ =
        unoptimizedController_.timeStamp_;

    for (size_t timeIndex = 0; timeIndex < PredictLength; ++timeIndex) {
      calculateControllerWorker(timeIndex, nominalPrimalData_, nominalDualData_,
                                unoptimizedController_);
    }

    // Since the controller for the last timestamp is invalid, if the last time
    // is not the event time, use the control policy of the second to last time
    // for the last time finalTimeIsNotAnEvent && there are at least two time
    // stamps
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
   * Calculate controller for the timeIndex by using primal and dual and write
   * the result back to dstController
   *
   * @param [in] timeIndex 当前时间索引。
   * @param [in] primalData 用于计算控制器的原始数据。
   * @param [in] dualData 用于计算控制器的对偶数据。
   * @param [out] dstController 输出的控制器（增益、偏置、deltaBias
   * 写入对应索引）。
   */
  void calculateControllerWorker(size_t timeIndex,
                                 const PrimalDataContainer_t& primalData,
                                 const DualDataContainer_t& dualData,
                                 LinearController_t& dstController) {
    const StateVector_t& nominalState =
        primalData.primalSolution.stateTrajectory_[timeIndex];
    const InputVector_t& nominalInput =
        primalData.primalSolution.inputTrajectory_[timeIndex];

    const Matrix<Scalar, UDim, UDim>& Qu =
        dualData.riccatiModificationTrajectory[timeIndex]
            .constraintNullProjector_;

    // feedback gains
    dstController.gainArray_[timeIndex] =
        Qu * projectedKmTrajectoryStock_[timeIndex];

    // bias input
    dstController.biasArray_[timeIndex] = nominalInput;
    dstController.biasArray_[timeIndex] -=
        dstController.gainArray_[timeIndex] * nominalState;
    dstController.deltaBiasArray_[timeIndex] =
        Qu * projectedLvTrajectoryStock_[timeIndex];
  }

  /** Based on the current LQ solution updates the optimized primal and dual
   * solutions. */
  void takePrimalDualStep(Scalar lqModelExpectedCost) {
    // update primal: run search strategy and find the optimal stepLength
    Scalar avgTimeStep = 0;
    SearchStrategySolutionRef_t solution(
        avgTimeStep, optimizedDualSolution_, optimizedPrimalSolution_,
        optimizedProblemMetrics_, performanceIndex_);
    const bool success = lineSearchStrategy_.run(
        {initTime_, finalTime_}, initState_, lqModelExpectedCost,
        unoptimizedController_, nominalDualData_.dualSolution, solution);

    if (success) {
      avgTimeStepFP_ = 0.9 * avgTimeStepFP_ + 0.1 * avgTimeStep;
    }

    // update dual
    if (success) {
      DualSolutionRef_t DualSolutionRef = optimizedDualSolution_;
      updateDualSolution(optimalControlProblem_, optimizedPrimalSolution_,
                         optimizedProblemMetrics_, DualSolutionRef);
      performanceIndex_ = computeRolloutPerformanceIndex(
          optimizedPrimalSolution_.timeTrajectory_, optimizedProblemMetrics_);
      performanceIndex_.merit = calculateRolloutMerit(performanceIndex_);
    }

    // if failed, use nominal and to keep the consistency of cached data, all
    // cache should be left untouched
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
    // result
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

    // Re-center around query state
    const StateVector_t xNominal = LinearInterpolation::interpolate(
        indexAlpha, primalSolution.stateTrajectory_);
    const StateVector_t deltaX = state - xNominal;
    const StateVector_t SmDeltaX = valueFunction.dfdxx * deltaX;
    valueFunction.f += deltaX.dot(0.5 * SmDeltaX + valueFunction.dfdx);
    valueFunction.dfdx += SmDeltaX;  // Adapt dfdx after f!

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
    valueFunction.dfdx += SmDeltaX;  // Adapt dfdx after f!

    return valueFunction;
  }

 public:
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
   * Calculates the merit function based on the performance index .
   *
   * @param [in] performanceIndex 性能指标（含 cost、等式/不等式拉格朗日等）。
   * @return 用于线搜索与收敛判据的 merit 值。
   */
  static Scalar calculateRolloutMerit(
      const PerformanceIndex_t& performanceIndex) {
    // cost
    Scalar merit = performanceIndex.cost;
    // state/state-input equality Lagrangian
    merit += performanceIndex.equalityLagrangian;
    // state/state-input inequality Lagrangian
    merit += performanceIndex.inequalityLagrangian;

    return merit;
  }

  /**
   * Computes cost, soft constraints and constraints values of each point in the
   * the primalSolution rollout.
   *
   * @param [in] problem 最优控制问题。
   * @param [in] primalSolution 原始解（轨迹）。
   * @param [in] dualSolution 对偶解。
   * @param [out] problemMetrics 各时刻代价、软约束与约束值。
   */
  static void computeRolloutMetrics(OptimalControlProblem_t& problem,
                                    const PrimalSolution_t& primalSolution,
                                    const DualSolution_t& dualSolution,
                                    ProblemMetrics_t& problemMetrics) {
    const TimeTrajectory_t& tTrajectory = primalSolution.timeTrajectory_;
    const StateTrajectory_t& xTrajectory = primalSolution.stateTrajectory_;
    const InputTrajectory_t& uTrajectory = primalSolution.inputTrajectory_;

    for (size_t k = 0; k < PredictLength; k++) {
      // intermediate time cost and constraints
      problemMetrics.intermediates[k] =
          approximator_.computeIntermediateMetrics(
              problem, tTrajectory[k], xTrajectory[k], uTrajectory[k],
              dualSolution.intermediates[k]);
    }

    // final time cost and constraints
    problemMetrics.final = approximator_.computeFinalMetrics(
        problem, tTrajectory.back(), xTrajectory.back(), dualSolution.final);
  }

  /**
   * Forward integrate the system dynamics with given controller. It uses the
   * given control policies and initial state, to integrate the system dynamics
   * in time period [initTime, finalTime].
   *
   * @param [in] rollout 前向积分用的 rollout 对象。
   * @param [in] initTime 初始时间。
   * @param [in] initState 初始状态。
   * @param [in] finalTime 终止时间。
   * @param [in,out] primalSolution 结果写入其时间/状态/输入轨迹；需已设置
   * controller 供 rollout 使用。
   * @return 平均时间步长。
   */
  static Scalar rolloutTrajectory(RolloutBase_t& rollout, Scalar initTime,
                                  const StateVector_t& initState,
                                  Scalar finalTime,
                                  PrimalSolution_t& primalSolution) {
    RolloutTrajectoryPointer_t rolloutTrajectortPtr(
        primalSolution.timeTrajectory_.data(),
        primalSolution.stateTrajectory_.data(),
        primalSolution.inputTrajectory_.data(), PredictLength + 1);
    rollout.run(initTime, initState, finalTime, &primalSolution.controller_,
                rolloutTrajectortPtr);
    // average time step
    return (finalTime - initTime) / static_cast<Scalar>(PredictLength);
  }

  /**
   * Calculates the PerformanceIndex associated to the given ProblemMetrics.
   *
   * @param [in] timeTrajectory rollout 的时间戳序列。
   * @param [in] problemMetrics 各时刻代价与约束指标。
   * @return 整条轨迹的 PerformanceIndex（梯形积分汇总）。
   */
  static PerformanceIndex_t computeRolloutPerformanceIndex(
      const TimeTrajectory_t& timeTrajectory,
      const ProblemMetrics_t& problemMetrics) {
    // Final
    PerformanceIndex_t finalperformanceIndex =
        toPerformanceIndex(problemMetrics.final);
    // Intermediates
    IntermediatePerformanceIndexTrajectory_t performanceIndexTrajectory;

    for (size_t i = 0; i < performanceIndexTrajectory.size(); ++i) {
      performanceIndexTrajectory[i] =
          toPerformanceIndex(problemMetrics.intermediates[i]);
    }

    // Intermediates
    return trapezoidalIntegration(timeTrajectory, performanceIndexTrajectory,
                                  finalperformanceIndex);
  }

  /**
   * @brief 按步长生成控制器：时间戳与增益与 unoptimizedController 相同，bias =
   * bias + stepLength * deltaBias。
   * @param [in] stepLength 线搜索步长。
   * @param [in] unoptimizedController 未优化控制器（含 deltaBiasArray）。
   * @param [out] controller 输出控制器。
   */
  static void incrementController(
      Scalar stepLength, const LinearController_t& unoptimizedController,
      LinearController_t& controller) {
    controller.timeStamp_ = unoptimizedController.timeStamp_;
    controller.gainArray_ = unoptimizedController.gainArray_;
    for (size_t k = 0; k < unoptimizedController.size(); k++) {
      controller.biasArray_[k] =
          unoptimizedController.biasArray_[k] +
          stepLength * unoptimizedController.deltaBiasArray_[k];
    }
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
   * Computes the integral of the squared (IS) norm of the controller update.
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
    // integrates using the trapezoidal approximation method
    return trapezoidalIntegration(controller.timeStamp_, biasArraySquaredNorm,
                                  0.0);
  }

 private:
  DDPSettings<Scalar> ddpSettings_{};

 public:
  OptimalControlProblem_t optimalControlProblem_;

  // roll out
  TimeTriggeredRollout_t rollout_;

  // initializer
  InitializerRollout_t initializerRollout_;

 private:
  // linear approximator
  inline static LinearQuadraticApproximator_t approximator_{};

  // time and state
  Scalar initTime_{0.0};
  Scalar finalTime_{0.0};
  Scalar lastFinalTime_{0.0};

  StateVector_t initState_;

  // nominal data
  PrimalDataContainer_t nominalPrimalData_;
  DualDataContainer_t nominalDualData_;

  // controller that is calculated directly from dual solution. It is
  // unoptimized because it haven't gone through searching.
  LinearController_t unoptimizedController_;

  LineSearchStrategy_t lineSearchStrategy_;

  // reference trajectory
  TimeTrajectory_t timeTrajectory_;
  StateTrajectory_t stateTrajectory_;
  InputTrajectory_t inputTrajectory_;

  // optimized data
  PrimalSolution_t optimizedPrimalSolution_;
  DualSolution_t optimizedDualSolution_;
  ProblemMetrics_t optimizedProblemMetrics_;

  // performance index
  PerformanceIndex_t performanceIndex_;
  PerformanceIndex_t performanceIndexLast_;

  // Discretizer
  EK2DynamicsDiscretizer_t discretizer_;

  // Discrete time riccati equation solver
  DiscreteTimeRiccatiEquations_t riccatiEquationsSolver_;
  std::array<KmMatrix_t, PredictLength + 1>
      projectedKmTrajectoryStock_;  // projected feedback
  std::array<LvVector_t, PredictLength + 1>
      projectedLvTrajectoryStock_;  // projected feedforward

  // forward pass and backward pass average time step
  Scalar avgTimeStepFP_ = 0.0;
  Scalar avgTimeStepBP_ = 0.0;

  size_t totalNumIterations_{0};
};