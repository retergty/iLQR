/**
 * @file LinearQuadraticApproximator.hpp
 * @brief 线性二次近似器：在名义轨迹上对最优控制问题做 LQ 近似，得到各节点
 * ModelData 与 Metrics。
 */
#pragma once

#include "LinearApproximation.hpp"
#include "Metrics.hpp"
#include "ModelData.hpp"
#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "QuadraticApproximation.hpp"
#include "iLQRDescriptorTraits.hpp"

/**
 * @brief 在给定名义轨迹与对偶解下，对 OCP 做 LQ 近似，填充中间/终端 ModelData
 * 与可选的 Metrics。
 * @tparam Descriptor iLQR 描述类型，提供状态/输入/时域和约束布局。
 */
template <typename Descriptor>
struct LinearQuadraticApproximator {
  using Traits = iLQRDescriptorTraits<Descriptor>;
  using Scalar = typename Traits::Scalar;
  using TranscriptionConfig = typename Traits::TranscriptionConfig;
  using Dims = typename Traits::Dims;
  using ConstraintConfig = typename Traits::ConstraintConfig;

  static constexpr int XDim = Traits::XDim;
  static constexpr int UDim = Traits::UDim;
  static constexpr std::size_t PredictLength = Traits::PredictLength;
  static constexpr int StateEqConstraintDim = Traits::StateEqDim;
  static constexpr int StateIneqConstraintDim = Traits::StateIneqDim;
  static constexpr int StateInputEqConstraintDim = Traits::StateInputEqDim;
  static constexpr int StateInputIneqConstraintDim = Traits::StateInputIneqDim;
  static constexpr int FinalStateEqConstraintDim = Traits::FinalStateEqDim;
  static constexpr int FinalStateIneqConstraintDim = Traits::FinalStateIneqDim;

  using OptimalControlProblem_t =
      OptimalControlProblem<Scalar, TranscriptionConfig, ConstraintConfig>;
  using StateVector_t = typename Traits::StateVector_t;
  using InputVector_t = typename Traits::InputVector_t;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;

  using IntermediateMultiplierCollection_t = MultiplierCollection<
      Scalar, typename Traits::IntermediateStageConstraintLayout_t>;
  using FinalMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           typename Traits::FinalStageConstraintLayout_t>;
  using IntermediateMetrics_t =
      Metrics<Scalar, Dims,
              typename Traits::IntermediateStageConstraintLayout_t>;
  using FinalMetrics_t =
      Metrics<Scalar, Dims, typename Traits::FinalStageConstraintLayout_t>;
  using TimeTrajectory_t = typename Traits::TimeTrajectory_t;
  using StateTrajectory_t = typename Traits::StateTrajectory_t;
  using InputTrajectory_t = typename Traits::InputTrajectory_t;
  using TargetTrajectories_t = typename Traits::TargetTrajectories_t;

  /**
   * 计算中间节点目标函数的二次近似，不计算动力学。
   *
   * 目标函数包括运行代价、状态代价以及增广拉格朗日项。该函数用于转录层
   * 自行生成连续/离散动力学近似时复用目标函数近似逻辑，避免先计算连续
   * 动力学线性化再被转录层覆盖。
   *
   * @param [in] problem 最优控制问题。
   * @param [in] targetTrajectories 参考轨迹。
   * @param [in] time 当前时间。
   * @param [in] state 当前状态。
   * @param [in] input 当前输入。
   * @param [in] multipliers 当前中间节点乘子。
   * @param [out] modelData 输出数据模型；仅写入 time 与 cost。
   */
  static void approximateIntermediateCost(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, const Scalar time,
      const StateVector_t& state, const InputVector_t& input,
      const IntermediateMultiplierCollection_t& multipliers,
      ModelData_t& modelData) {
    modelData.time = time;

    // 代价。
    modelData.cost =
        approximateCost(problem, targetTrajectories, time, state, input);

    // 拉格朗日项。
    if constexpr (StateEqConstraintDim != 0) {
      ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> approx =
          problem.stateEqualityLagrangian.getQuadraticApproximation(
              time, state, multipliers.stateEq);
      modelData.cost.f += approx.f;
      modelData.cost.dfdx += approx.dfdx;
      modelData.cost.dfdxx += approx.dfdxx;
    }
    if constexpr (StateIneqConstraintDim != 0) {
      ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> approx =
          problem.stateInequalityLagrangian.getQuadraticApproximation(
              time, state, multipliers.stateIneq);
      modelData.cost.f += approx.f;
      modelData.cost.dfdx += approx.dfdx;
      modelData.cost.dfdxx += approx.dfdxx;
    }
    if constexpr (StateInputEqConstraintDim != 0) {
      modelData.cost += problem.equalityLagrangian.getQuadraticApproximation(
          time, state, input, multipliers.stateInputEq);
    }
    if constexpr (StateInputIneqConstraintDim != 0) {
      modelData.cost += problem.inequalityLagrangian.getQuadraticApproximation(
          time, state, input, multipliers.stateInputIneq);
    }
  }

  /**
   * 计算中间节点目标函数的二次近似，不计算动力学。
   *
   * @return 仅包含 time 与 cost 的模型数据。
   */
  static inline ModelData_t approximateIntermediateCost(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, const Scalar time,
      const StateVector_t& state, const InputVector_t& input,
      const IntermediateMultiplierCollection_t& multipliers) {
    ModelData_t md;
    approximateIntermediateCost(problem, targetTrajectories, time, state, input,
                                multipliers, md);
    return md;
  }

  /**
   * @brief 计算终端节点的目标函数与终端增广拉格朗日二次近似。
   *
   * @param [in] problem: 最优控制问题。
   * @param [in] time: 当前时间。
   * @param [in] state: 当前状态。
   * @param [in] multipliers: 当前与等式和不等式拉格朗日项关联的乘子。
   * @param [out] modelData: 输出数据模型。
   */
  static void approximateFinalLQ(const OptimalControlProblem_t& problem,
                                 const TargetTrajectories_t& targetTrajectories,
                                 const Scalar time, const StateVector_t& state,
                                 const FinalMultiplierCollection_t& multipliers,
                                 ModelData_t& modelData) {
    modelData.time = time;

    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> finalDynamics;
    finalDynamics.setZero();

    // 终端节点没有后续动力学，填零以保持 ModelData 结构完整。
    modelData.dynamics = finalDynamics;

    // 终端代价
    modelData.cost =
        approximateFinalCost(problem, targetTrajectories, time, state);

    // 拉格朗日项
    if constexpr (FinalStateEqConstraintDim != 0) {
      auto approx = problem.finalEqualityLagrangian.getQuadraticApproximation(
          time, state, multipliers.stateEq);
      modelData.cost.f += approx.f;
      modelData.cost.dfdx += approx.dfdx;
      modelData.cost.dfdxx += approx.dfdxx;
    }
    if constexpr (FinalStateIneqConstraintDim != 0) {
      auto approx = problem.finalInequalityLagrangian.getQuadraticApproximation(
          time, state, multipliers.stateIneq);
      modelData.cost.f += approx.f;
      modelData.cost.dfdx += approx.dfdx;
      modelData.cost.dfdxx += approx.dfdxx;
    }
  }

  /**
   * @brief 返回终端节点的目标函数与终端增广拉格朗日二次近似。
   *
   * @param [in] problem: 最优控制问题。
   * @param [in] time: 当前时间。
   * @param [in] state: 当前状态。
   * @param [in] multipliers: 当前与等式和不等式拉格朗日项关联的乘子。
   * @return 输出数据模型。
   */
  static inline ModelData_t approximateFinalLQ(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, const Scalar time,
      const StateVector_t& state,
      const FinalMultiplierCollection_t& multipliers) {
    ModelData_t md;
    approximateFinalLQ(problem, targetTrajectories, time, state, multipliers,
                       md);
    return md;
  }

  /**
   * 计算中间时刻总代价（即 cost + softConstraints）。
   * 假定预计算请求已经发出。
   */
  static Scalar computeCost(const OptimalControlProblem_t& problem,
                            const TargetTrajectories_t& targetTrajectories,
                            const Scalar time, const StateVector_t& state,
                            const InputVector_t& input) {
    const TimeTrajectory_t& targetTimeTrajectories =
        targetTrajectories.timeTrajectory;
    const StateTrajectory_t& targetStateTrajectories =
        targetTrajectories.stateTrajectory;
    const InputTrajectory_t& targetInputTrajectories =
        targetTrajectories.inputTrajectory;

    // 计算并累加所有代价。
    Scalar cost =
        problem.cost.getValue(time, state, input, targetTimeTrajectories,
                              targetStateTrajectories, targetInputTrajectories);
    cost += problem.stateCost.getValue(time, state, targetTimeTrajectories,
                                       targetStateTrajectories);

    return cost;
  }

  /**
   * 计算中间时刻总代价（即
   * cost + softConstraints）的二次近似。假定预计算请求
   * 已经发出。
   */
  static ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  approximateCost(const OptimalControlProblem_t& problem,
                  const TargetTrajectories_t& targetTrajectories,
                  const Scalar time, const StateVector_t& state,
                  const InputVector_t& input) {
    const TimeTrajectory_t& targetTimeTrajectories =
        targetTrajectories.timeTrajectory;
    const StateTrajectory_t& targetStateTrajectories =
        targetTrajectories.stateTrajectory;
    const InputTrajectory_t& targetInputTrajectories =
        targetTrajectories.inputTrajectory;

    // 获取状态-输入代价近似。
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> cost =
        problem.cost.getQuadraticApproximation(
            time, state, input, targetTimeTrajectories, targetStateTrajectories,
            targetInputTrajectories);

    // 获取仅状态代价近似。
    cost += problem.stateCost.getQuadraticApproximation(
        time, state, targetTimeTrajectories, targetStateTrajectories);

    return cost;
  }

  /**
   * 计算终端总代价（即 cost + softConstraints）。假定
   * 预计算请求已经发出。
   */
  static Scalar computeFinalCost(const OptimalControlProblem_t& problem,
                                 const TargetTrajectories_t& targetTrajectories,
                                 const Scalar time,
                                 const StateVector_t& state) {
    const TimeTrajectory_t& targetTimeTrajectories =
        targetTrajectories.timeTrajectory;
    const StateTrajectory_t& targetStateTrajectories =
        targetTrajectories.stateTrajectory;

    Scalar cost = problem.finalCost.getValue(
        time, state, targetTimeTrajectories, targetStateTrajectories);

    return cost;
  }

  /**
   * 计算终端总代价（即 cost +
   * softConstraints）的二次近似。假定预计算请求已经
   * 发出。
   */
  static ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  approximateFinalCost(const OptimalControlProblem_t& problem,
                       const TargetTrajectories_t& targetTrajectories,
                       const Scalar time, const StateVector_t& state) {
    const TimeTrajectory_t& targetTimeTrajectories =
        targetTrajectories.timeTrajectory;
    const StateTrajectory_t& targetStateTrajectories =
        targetTrajectories.stateTrajectory;

    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> cost =
        problem.finalCost.getQuadraticApproximation(
            time, state, targetTimeTrajectories, targetStateTrajectories);

    return cost;
  }

  /**
   * 计算中间时刻 Metrics（即 cost、softConstraints 和
   * 约束）。
   *
   * @note 假定预计算请求已经发出。
   * problem.preComputationPtr->request(Request::Cost + Request::Constraint +
   * Request::SoftConstraint, t, x, u)
   *
   * @param [in] problem: 最优控制问题。
   * @param [in] time: 当前时间。
   * @param [in] state: 当前状态。
   * @param [in] input: 当前输入。
   * @param [in] dynamicsViolation: 动力学违反量，取决于
   * 转录方法。
   * @return 输出 Metrics。
   */
  static IntermediateMetrics_t computeIntermediateMetrics(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, const Scalar time,
      const StateVector_t& state, const InputVector_t& input) {
    IntermediateMetrics_t metrics;

    // 代价
    metrics.cost = computeCost(problem, targetTrajectories, time, state, input);

    return metrics;
  }

  /**
   * 计算中间时刻 Metrics（即 cost、softConstraints 和
   * 约束）。
   *
   * @note 假定预计算请求已经发出。
   * problem.preComputationPtr->request(Request::Cost + Request::Constraint +
   * Request::SoftConstraint, t, x, u)
   *
   * @param [in] problem: 最优控制问题。
   * @param [in] time: 当前时间。
   * @param [in] state: 当前状态。
   * @param [in] input: 当前输入。
   * @param [in] multipliers: 当前与等式和不等式拉格朗日项关联的乘子。
   * @param [in] dynamicsViolation: 动力学违反量，取决于
   * 转录方法。
   * @return 输出 Metrics。
   */
  static IntermediateMetrics_t computeIntermediateMetrics(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, const Scalar time,
      const StateVector_t& state, const InputVector_t& input,
      const IntermediateMultiplierCollection_t& multipliers) {
    // 代价、动力学违反、等式约束和不等式约束。
    IntermediateMetrics_t metrics = computeIntermediateMetrics(
        problem, targetTrajectories, time, state, input);

    // 等式拉格朗日项
    if constexpr (StateEqConstraintDim != 0) {
      metrics.stateEqLagrangian = problem.stateEqualityLagrangian.getValue(
          time, state, multipliers.stateEq);
    }

    if constexpr (StateInputEqConstraintDim != 0) {
      metrics.stateInputEqLagrangian = problem.equalityLagrangian.getValue(
          time, state, input, multipliers.stateInputEq);
    }

    // 不等式拉格朗日项
    if constexpr (StateIneqConstraintDim != 0) {
      metrics.stateIneqLagrangian = problem.stateInequalityLagrangian.getValue(
          time, state, multipliers.stateIneq);
    }

    if constexpr (StateInputIneqConstraintDim != 0) {
      metrics.stateInputIneqLagrangian = problem.inequalityLagrangian.getValue(
          time, state, input, multipliers.stateInputIneq);
    }

    return metrics;
  }

  /**
   * 计算终端时刻 Metrics（即 cost、softConstraints 和
   * 约束）。
   *
   * @note 假定预计算请求已经发出。
   * problem.preComputationPtr->requestFinal(Request::Cost + Request::Constraint
   * + Request::SoftConstraint, t, x)
   *
   * @param [in] problem: 最优控制问题。
   * @param [in] time: 当前时间。
   * @param [in] state: 当前状态。
   * @return 输出 Metrics。
   */
  static FinalMetrics_t computeFinalMetrics(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, const Scalar time,
      const StateVector_t& state) {
    FinalMetrics_t metrics;

    // 代价
    metrics.cost = computeFinalCost(problem, targetTrajectories, time, state);

    return metrics;
  }

  /**
   * 计算终端时刻 Metrics（即 cost、softConstraints 和
   * 约束）。
   *
   * @note 假定预计算请求已经发出。
   * problem.preComputationPtr->requestFinal(Request::Cost + Request::Constraint
   * + Request::SoftConstraint, t, x)
   *
   * @param [in] problem: 最优控制问题。
   * @param [in] time: 当前时间。
   * @param [in] state: 当前状态。
   * @param [in] multipliers: 当前与等式和不等式拉格朗日项关联的乘子。
   * @return 输出 Metrics。
   */
  static FinalMetrics_t computeFinalMetrics(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, const Scalar time,
      const StateVector_t& state,
      const FinalMultiplierCollection_t& multipliers) {
    // 代价、等式约束和不等式约束。
    FinalMetrics_t metrics =
        computeFinalMetrics(problem, targetTrajectories, time, state);

    if constexpr (FinalStateEqConstraintDim != 0) {
      // 等式拉格朗日项
      metrics.stateEqLagrangian = problem.finalEqualityLagrangian.getValue(
          time, state, multipliers.stateEq);
    }

    // 不等式拉格朗日项
    if constexpr (FinalStateIneqConstraintDim != 0) {
      metrics.stateIneqLagrangian = problem.finalInequalityLagrangian.getValue(
          time, state, multipliers.stateIneq);
    }

    return metrics;
  }
};
