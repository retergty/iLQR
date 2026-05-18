/**
 * @file OptimalControlProblem.hpp
 * @brief 最优控制问题定义：代价、拉格朗日、参考轨迹与动力学指针。
 */
#pragma once

#include "ControlledSystemBase.hpp"
#include "Controller.hpp"
#include "Cost.hpp"
#include "CostCollection.hpp"
#include "Multiplier.hpp"
#include "StateAugmentedLagrangianCollection.hpp"
#include "StateInputAugmentedLagrangianCollection.hpp"
#include "SystemDynamicsBase.hpp"

/**
 * @brief
 * 最优控制问题：中间/终端代价、等式/不等式增广拉格朗日、参考轨迹与系统动力学。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 * @tparam PredictLength 预测步数。
 * @tparam StateEqLagrangianContrainNumbers 等 各约束维度。
 */
template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianContrainNumbers,
          int StateInputEqLagrangianContrainNumbers,
          int StateIneqLagrangianContrainNumbers,
          int StateInputIneqLagrangianContrainNumbers,
          int FinalStateEqLagrangianContrainNumbers,
          int FinalStateIneqFinalLagrangianContrainNumbers>
struct OptimalControlProblem {
  /** @brief 默认构造。 */
  OptimalControlProblem() = default;

  /** @brief 析构函数。 */
  ~OptimalControlProblem() = default;

  /** @brief 禁止拷贝构造。 */
  OptimalControlProblem(const OptimalControlProblem& other) = delete;

  /** @brief 禁止拷贝赋值。 */
  OptimalControlProblem& operator=(const OptimalControlProblem& rhs) = delete;

  // /** Move constructor */
  // OptimalControlProblem(OptimalControlProblem &&other) noexcept = default;

  // /** Move assignment */
  // OptimalControlProblem &operator=(OptimalControlProblem &&rhs) noexcept =
  // default;

  /** @brief 中间代价（状态-输入）。 */
  StateInputCostCollection<Scalar, XDim, UDim, PredictLength + 1> cost;
  /** @brief 仅状态中间代价。 */
  StateCostCollection<Scalar, XDim, PredictLength + 1> stateCost;

  /** @brief 终端代价。 */
  StateCostCollection<Scalar, XDim, PredictLength + 1> finalCost;

  /** @brief 参考时间轨迹。 */
  std::array<Scalar, PredictLength + 1> timeTrajectory;
  /** @brief 参考状态轨迹。 */
  std::array<Vector<Scalar, XDim>, PredictLength + 1> stateTrajectory;
  /** @brief 参考输入轨迹。 */
  std::array<Vector<Scalar, UDim>, PredictLength + 1> inputTrajectory;

  /** @brief 状态-输入等式约束增广拉格朗日。 */
  StateInputAugmentedLagrangianCollection<Scalar, XDim, UDim,
                                          StateInputEqLagrangianContrainNumbers>
      equalityLagrangian;
  /** @brief 仅状态等式约束增广拉格朗日。 */
  StateAugmentedLagrangianCollection<Scalar, XDim,
                                     StateEqLagrangianContrainNumbers>
      stateEqualityLagrangian;
  /** @brief 状态-输入不等式约束增广拉格朗日。 */
  StateInputAugmentedLagrangianCollection<
      Scalar, XDim, UDim, StateInputIneqLagrangianContrainNumbers>
      inequalityLagrangian;
  /** @brief 仅状态不等式约束增广拉格朗日。 */
  StateAugmentedLagrangianCollection<Scalar, XDim,
                                     StateIneqLagrangianContrainNumbers>
      stateInequalityLagrangian;

  /** @brief 终端等式约束增广拉格朗日。 */
  StateAugmentedLagrangianCollection<Scalar, XDim,
                                     FinalStateEqLagrangianContrainNumbers>
      finalEqualityLagrangian;
  /** @brief 终端不等式约束增广拉格朗日。 */
  StateAugmentedLagrangianCollection<
      Scalar, XDim, FinalStateIneqFinalLagrangianContrainNumbers>
      finalInequalityLagrangian;

  /** @brief 系统动力学指针。 */
  SystemDynamicsBase<Scalar, XDim, UDim>* dynamicsPtr;
};
