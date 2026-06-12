/**
 * @file OptimalControlProblem.hpp
 * @brief 最优控制问题定义：代价、拉格朗日、参考轨迹与动力学指针。
 */
#pragma once

#include "AugmentedLagrangian/StateAugmentedLagrangianCollection.hpp"
#include "AugmentedLagrangian/StateInputAugmentedLagrangianCollection.hpp"
#include "Cost/CostCollection.hpp"
#include "Dynamics/DynamicsModeTraits.hpp"

/**
 * @brief
 * 最优控制问题：中间/终端代价、等式/不等式增广拉格朗日、参考轨迹与系统动力学。
 * @tparam Scalar 标量类型。
 * @tparam Transcription 轨迹配置，提供 XDim/UDim/PredictLength。
 * @tparam ConstraintConfig 约束配置，提供各类约束的 term 分组布局与总维度。
 */
template <typename Scalar, typename Transcription, typename ConstraintConfig>
struct OptimalControlProblem {
  static constexpr int XDim = Transcription::XDim;
  static constexpr int UDim = Transcription::UDim;
  static constexpr std::size_t PredictLength = Transcription::PredictLength;

  using DynamicsMode = typename Transcription::DynamicsMode;
  using DynamicsBase_t =
      typename DynamicsModeTraits<Scalar, XDim, UDim,
                                  DynamicsMode>::DynamicsBase_t;

  using StateEqLayout = typename ConstraintConfig::StateEqLayout;
  using StateIneqLayout = typename ConstraintConfig::StateIneqLayout;
  using StateInputEqLayout = typename ConstraintConfig::StateInputEqLayout;
  using StateInputIneqLayout = typename ConstraintConfig::StateInputIneqLayout;
  using FinalStateEqLayout = typename ConstraintConfig::FinalStateEqLayout;
  using FinalStateIneqLayout = typename ConstraintConfig::FinalStateIneqLayout;

  static constexpr int StateEqConstraintDim = StateEqLayout::TotalDim;
  static constexpr int StateIneqConstraintDim = StateIneqLayout::TotalDim;
  static constexpr int StateInputEqConstraintDim = StateInputEqLayout::TotalDim;
  static constexpr int StateInputIneqConstraintDim =
      StateInputIneqLayout::TotalDim;
  static constexpr int FinalStateEqConstraintDim = FinalStateEqLayout::TotalDim;
  static constexpr int FinalStateIneqConstraintDim =
      FinalStateIneqLayout::TotalDim;

  /** @brief 默认构造。 */
  OptimalControlProblem() = default;

  /** @brief 析构函数。 */
  ~OptimalControlProblem() = default;

  /** @brief 禁止拷贝构造。 */
  OptimalControlProblem(const OptimalControlProblem& other) = delete;

  /** @brief 禁止拷贝赋值。 */
  OptimalControlProblem& operator=(const OptimalControlProblem& rhs) = delete;

  /** @brief 中间代价（状态-输入）。 */
  StateInputCostCollection<Scalar, XDim, UDim, PredictLength + 1> cost;
  /** @brief 仅状态中间代价。 */
  StateCostCollection<Scalar, XDim, UDim, PredictLength + 1> stateCost;

  /** @brief 终端代价。 */
  StateCostCollection<Scalar, XDim, UDim, PredictLength + 1> finalCost;

  /** @brief 状态-输入等式约束增广拉格朗日。 */
  StateInputAugmentedLagrangianCollection<Scalar, XDim, UDim,
                                          StateInputEqLayout>
      equalityLagrangian;
  /** @brief 仅状态等式约束增广拉格朗日。 */
  StateAugmentedLagrangianCollection<Scalar, XDim, UDim, StateEqLayout>
      stateEqualityLagrangian;
  /** @brief 状态-输入不等式约束增广拉格朗日。 */
  StateInputAugmentedLagrangianCollection<Scalar, XDim, UDim,
                                          StateInputIneqLayout>
      inequalityLagrangian;
  /** @brief 仅状态不等式约束增广拉格朗日。 */
  StateAugmentedLagrangianCollection<Scalar, XDim, UDim, StateIneqLayout>
      stateInequalityLagrangian;

  /** @brief 终端等式约束增广拉格朗日。 */
  StateAugmentedLagrangianCollection<Scalar, XDim, UDim, FinalStateEqLayout>
      finalEqualityLagrangian;
  /** @brief 终端不等式约束增广拉格朗日。 */
  StateAugmentedLagrangianCollection<Scalar, XDim, UDim, FinalStateIneqLayout>
      finalInequalityLagrangian;

  /** @brief 系统动力学指针，类型由 Transcription::DynamicsMode 决定。 */
  DynamicsBase_t* dynamicsPtr{nullptr};
};
