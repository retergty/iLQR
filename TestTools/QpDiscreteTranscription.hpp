/**
 * @file QpDiscreteTranscription.hpp
 * @brief 将名义轨迹附近的 OCP 转录为稠密 QP 测试工具使用的离散 LQ 问题。
 *
 * 中间阶段的动力学与目标函数近似复用主项目 Transcription 层，因此连续系统和
 * 离散系统会分别使用与 iLQR 主求解器一致的转录语义。该文件仍然使用动态尺寸
 * QP 数据结构，仅作为测试和对照工具，不属于核心实时路径。
 */
#pragma once

#include <array>
#include <stdexcept>
#include <vector>

#include "LinearQuadraticApproximator.hpp"
#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "QpSolverTypes.hpp"
#include "QpTrajectories.hpp"
#include "Reference.hpp"
#include "StateInputConstraint.hpp"
#include "Transcription/TranscriptionTraits.hpp"
#include "iLQRDescriptor.hpp"

namespace qp_solver {

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpContinuousDynamicsTranscriptionConfig_t =
    TranscriptionConfig<Dimensions<XDim, UDim>, Horizon<PredictLength>,
                        ContinuousDynamics>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpDiscreteDynamicsTranscriptionConfig_t =
    TranscriptionConfig<Dimensions<XDim, UDim>, Horizon<PredictLength>,
                        DiscreteDynamics>;

using QpConstraintConfig_t = ConstraintConfig<>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpContinuousDynamicsDescriptor_t =
    iLQRDescriptor<Scalar,
                   QpContinuousDynamicsTranscriptionConfig_t<Scalar, XDim, UDim,
                                                             PredictLength>,
                   QpConstraintConfig_t>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpDiscreteDynamicsDescriptor_t = iLQRDescriptor<
    Scalar,
    QpDiscreteDynamicsTranscriptionConfig_t<Scalar, XDim, UDim, PredictLength>,
    QpConstraintConfig_t>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpContinuousDynamicsOptimalControlProblem_t =
    typename LinearQuadraticApproximator<QpContinuousDynamicsDescriptor_t<
        Scalar, XDim, UDim, PredictLength>>::OptimalControlProblem_t;

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpDiscreteDynamicsOptimalControlProblem_t =
    typename LinearQuadraticApproximator<QpDiscreteDynamicsDescriptor_t<
        Scalar, XDim, UDim, PredictLength>>::OptimalControlProblem_t;

template <typename Scalar>
using QpIntermediateMultiplierCollection_t = MultiplierCollection<
    Scalar, IntermediateStageConstraintLayout<QpConstraintConfig_t>>;

template <typename Scalar>
using QpFinalMultiplierCollection_t =
    MultiplierCollection<Scalar,
                         FinalStageConstraintLayout<QpConstraintConfig_t>>;

template <typename ConstraintConfig>
constexpr bool isUnconstrainedQpConfig =
    ConstraintConfig::StateEq == 0 && ConstraintConfig::StateIneq == 0 &&
    ConstraintConfig::StateInputEq == 0 &&
    ConstraintConfig::StateInputIneq == 0 &&
    ConstraintConfig::FinalStateEq == 0 &&
    ConstraintConfig::FinalStateIneq == 0;

template <typename Scalar, int XDim, int UDim>
typename LinearQuadraticStage<Scalar, XDim, UDim>::ConstraintApproximation_t
makeZeroConstraints() {
  typename LinearQuadraticStage<Scalar, XDim, UDim>::ConstraintApproximation_t
      constraints;
  constraints.f.resize(0);
  constraints.dfdx.resize(0, XDim);
  constraints.dfdu.resize(0, UDim);
  return constraints;
}

template <typename Scalar, int XDim, int UDim>
typename LinearQuadraticStage<Scalar, XDim, UDim>::DynamicsApproximation_t
makeZeroDynamics() {
  typename LinearQuadraticStage<Scalar, XDim, UDim>::DynamicsApproximation_t
      dynamics;
  dynamics.setZero();
  return dynamics;
}

template <typename Scalar, typename Transcription>
TargetTrajectories<Scalar, Transcription> toTargetTrajectories(
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& trajectory) {
  TargetTrajectories<Scalar, Transcription> targetTrajectory;
  targetTrajectory.timeTrajectory = trajectory.timeTrajectory;
  targetTrajectory.stateTrajectory = trajectory.stateTrajectory;
  for (size_t k = 0; k < Transcription::PredictLength; ++k) {
    targetTrajectory.inputTrajectory[k] = trajectory.inputTrajectory[k];
  }
  if constexpr (Transcription::PredictLength > 0) {
    targetTrajectory.inputTrajectory[Transcription::PredictLength] =
        trajectory.inputTrajectory[Transcription::PredictLength - 1];
  }
  return targetTrajectory;
}

template <typename Scalar, int XDim, int UDim, int CDim>
typename LinearQuadraticStage<Scalar, XDim, UDim>::ConstraintApproximation_t
makeStateInputConstraints(
    const StateInputConstraint<Scalar, XDim, UDim, CDim>& constraint,
    TrajectoryRef<Scalar, XDim, UDim> start) {
  const auto linearConstraint =
      constraint.getLinearApproximation(start.t, start.x, start.u);
  typename LinearQuadraticStage<Scalar, XDim, UDim>::ConstraintApproximation_t
      constraints;
  constraints.f = linearConstraint.f;
  constraints.dfdx = linearConstraint.dfdx;
  constraints.dfdu = linearConstraint.dfdu;
  return constraints;
}

template <typename Scalar, typename Transcription, typename ConstraintConfig>
LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>
approximateStage(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const TargetTrajectories<Scalar, Transcription>& targetTrajectory,
    TrajectoryRef<Scalar, Transcription::XDim, Transcription::UDim> start,
    StateTrajectoryRef<Scalar, Transcription::XDim> end,
    const MultiplierCollection<
        Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>&
        multipliers) {
  static_assert(isUnconstrainedQpConfig<ConstraintConfig>,
                "Dense QP transcription only supports unconstrained optimal "
                "control problems. Use ConstraintConfig<> and keep augmented "
                "Lagrangian dimensions at zero.");

  constexpr int XDim = Transcription::XDim;
  constexpr int UDim = Transcription::UDim;

  using Descriptor_t = iLQRDescriptor<Scalar, Transcription, ConstraintConfig>;
  using TranscriptionImpl_t = ::Transcription_t<Descriptor_t>;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;
  using LqStage_t = LinearQuadraticStage<Scalar, XDim, UDim>;

  LqStage_t lqStage;
  const Scalar dt = end.t - start.t;
  TranscriptionImpl_t transcription;
  ModelData_t modelData;
  transcription.approximateIntermediateLQ(optimalControlProblem,
                                          targetTrajectory, start.t, start.x,
                                          start.u, dt, multipliers, modelData);

  lqStage.cost = modelData.cost;
  lqStage.dynamics = modelData.dynamics;

  // 增广拉格朗日约束惩罚已经由主项目 Transcription/Approximator
  // 折叠进目标函数近似；这里不再生成额外软约束。
  lqStage.constraints = makeZeroConstraints<Scalar, XDim, UDim>();

  return lqStage;
}

template <typename Scalar, typename Transcription, typename ConstraintConfig,
          int CDim>
LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>
approximateStage(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const TargetTrajectories<Scalar, Transcription>& targetTrajectory,
    TrajectoryRef<Scalar, Transcription::XDim, Transcription::UDim> start,
    StateTrajectoryRef<Scalar, Transcription::XDim> end,
    const MultiplierCollection<
        Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>&
        multipliers,
    const StateInputConstraint<Scalar, Transcription::XDim, Transcription::UDim,
                               CDim>& constraint) {
  constexpr int XDim = Transcription::XDim;
  constexpr int UDim = Transcription::UDim;

  auto lqStage = approximateStage(optimalControlProblem, targetTrajectory,
                                  start, end, multipliers);
  lqStage.constraints =
      makeStateInputConstraints<Scalar, XDim, UDim>(constraint, start);
  return lqStage;
}

template <typename Scalar, typename Transcription, typename ConstraintConfig>
std::vector<
    LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& nominalTrajectory,
    const std::array<
        MultiplierCollection<
            Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
        Transcription::PredictLength>& intermediateMultipliers,
    const MultiplierCollection<Scalar,
                               FinalStageConstraintLayout<ConstraintConfig>>&
        finalMultipliers) {
  const auto targetTrajectory =
      toTargetTrajectories<Scalar, Transcription>(nominalTrajectory);
  return getLinearQuadraticApproximation(
      optimalControlProblem, targetTrajectory, nominalTrajectory,
      intermediateMultipliers, finalMultipliers);
}

template <typename Scalar, typename Transcription, typename ConstraintConfig>
std::vector<
    LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const TargetTrajectories<Scalar, Transcription>& targetTrajectory,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& nominalTrajectory,
    const std::array<
        MultiplierCollection<
            Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
        Transcription::PredictLength>& intermediateMultipliers,
    const MultiplierCollection<Scalar,
                               FinalStageConstraintLayout<ConstraintConfig>>&
        finalMultipliers) {
  static_assert(isUnconstrainedQpConfig<ConstraintConfig>,
                "Dense QP transcription only supports unconstrained optimal "
                "control problems. Use ConstraintConfig<> and keep augmented "
                "Lagrangian dimensions at zero.");

  constexpr int XDim = Transcription::XDim;
  constexpr int UDim = Transcription::UDim;
  constexpr size_t PredictLength = Transcription::PredictLength;

  if (optimalControlProblem.dynamicsPtr == nullptr) {
    throw std::runtime_error(
        "[qp_solver::getLinearQuadraticApproximation] "
        "dynamicsPtr should not be null.");
  }

  using Descriptor_t = iLQRDescriptor<Scalar, Transcription, ConstraintConfig>;
  using TranscriptionImpl_t = ::Transcription_t<Descriptor_t>;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;

  const auto& t = nominalTrajectory.timeTrajectory;
  const auto& x = nominalTrajectory.stateTrajectory;
  const auto& u = nominalTrajectory.inputTrajectory;
  constexpr size_t N = PredictLength;

  // 包含 N+1 个元素的 LinearQuadraticProblem；终端阶段 lqp[N].dynamics 会被
  // QP 构造过程忽略。
  std::vector<LinearQuadraticStage<Scalar, XDim, UDim>> lqp;
  lqp.reserve(N + 1);
  for (size_t k = 0; k < N; ++k) {  // 中间阶段。
    lqp.emplace_back(approximateStage(optimalControlProblem, targetTrajectory,
                                      {t[k], x[k], u[k]}, {t[k + 1], x[k + 1]},
                                      intermediateMultipliers[k]));
  }

  TranscriptionImpl_t transcription;
  ModelData_t modelData;
  transcription.approximateFinalLQ(optimalControlProblem, targetTrajectory,
                                   t[N], x[N], finalMultipliers, modelData);
  lqp.emplace_back(std::move(modelData.cost),
                   makeZeroDynamics<Scalar, XDim, UDim>(),
                   makeZeroConstraints<Scalar, XDim, UDim>());

  return lqp;
}

template <typename Scalar, typename Transcription, typename ConstraintConfig,
          int CDim>
std::vector<
    LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const TargetTrajectories<Scalar, Transcription>& targetTrajectory,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& nominalTrajectory,
    const std::array<
        MultiplierCollection<
            Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
        Transcription::PredictLength>& intermediateMultipliers,
    const MultiplierCollection<
        Scalar, FinalStageConstraintLayout<ConstraintConfig>>& finalMultipliers,
    const StateInputConstraint<Scalar, Transcription::XDim, Transcription::UDim,
                               CDim>& constraint) {
  static_assert(isUnconstrainedQpConfig<ConstraintConfig>,
                "Dense QP transcription only supports unconstrained optimal "
                "control problems. Use ConstraintConfig<> and keep augmented "
                "Lagrangian dimensions at zero.");

  constexpr int XDim = Transcription::XDim;
  constexpr int UDim = Transcription::UDim;
  constexpr size_t PredictLength = Transcription::PredictLength;

  if (optimalControlProblem.dynamicsPtr == nullptr) {
    throw std::runtime_error(
        "[qp_solver::getLinearQuadraticApproximation] "
        "dynamicsPtr should not be null.");
  }

  using Descriptor_t = iLQRDescriptor<Scalar, Transcription, ConstraintConfig>;
  using TranscriptionImpl_t = ::Transcription_t<Descriptor_t>;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;

  const auto& t = nominalTrajectory.timeTrajectory;
  const auto& x = nominalTrajectory.stateTrajectory;
  const auto& u = nominalTrajectory.inputTrajectory;
  constexpr size_t N = PredictLength;

  std::vector<LinearQuadraticStage<Scalar, XDim, UDim>> lqp;
  lqp.reserve(N + 1);
  for (size_t k = 0; k < N; ++k) {
    lqp.emplace_back(approximateStage(optimalControlProblem, targetTrajectory,
                                      {t[k], x[k], u[k]}, {t[k + 1], x[k + 1]},
                                      intermediateMultipliers[k], constraint));
  }

  TranscriptionImpl_t transcription;
  ModelData_t modelData;
  transcription.approximateFinalLQ(optimalControlProblem, targetTrajectory,
                                   t[N], x[N], finalMultipliers, modelData);
  lqp.emplace_back(std::move(modelData.cost),
                   makeZeroDynamics<Scalar, XDim, UDim>(),
                   makeZeroConstraints<Scalar, XDim, UDim>());

  return lqp;
}

template <typename Scalar, typename Transcription, typename ConstraintConfig>
std::vector<
    LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const TargetTrajectories<Scalar, Transcription>& targetTrajectory,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>&
        nominalTrajectory) {
  std::array<MultiplierCollection<
                 Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
             Transcription::PredictLength>
      intermediateMultipliers{};
  MultiplierCollection<Scalar, FinalStageConstraintLayout<ConstraintConfig>>
      finalMultipliers{};
  return getLinearQuadraticApproximation(
      optimalControlProblem, targetTrajectory, nominalTrajectory,
      intermediateMultipliers, finalMultipliers);
}

template <typename Scalar, typename Transcription, typename ConstraintConfig,
          int CDim>
std::vector<
    LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const TargetTrajectories<Scalar, Transcription>& targetTrajectory,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& nominalTrajectory,
    const StateInputConstraint<Scalar, Transcription::XDim, Transcription::UDim,
                               CDim>& constraint) {
  std::array<MultiplierCollection<
                 Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
             Transcription::PredictLength>
      intermediateMultipliers{};
  MultiplierCollection<Scalar, FinalStageConstraintLayout<ConstraintConfig>>
      finalMultipliers{};
  return getLinearQuadraticApproximation(
      optimalControlProblem, targetTrajectory, nominalTrajectory,
      intermediateMultipliers, finalMultipliers, constraint);
}

template <typename Scalar, typename Transcription, typename ConstraintConfig>
std::vector<
    LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>&
        nominalTrajectory) {
  const auto targetTrajectory =
      toTargetTrajectories<Scalar, Transcription>(nominalTrajectory);
  std::array<MultiplierCollection<
                 Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
             Transcription::PredictLength>
      intermediateMultipliers{};
  MultiplierCollection<Scalar, FinalStageConstraintLayout<ConstraintConfig>>
      finalMultipliers{};
  return getLinearQuadraticApproximation(
      optimalControlProblem, targetTrajectory, nominalTrajectory,
      intermediateMultipliers, finalMultipliers);
}

template <typename Scalar, typename Transcription, typename ConstraintConfig,
          int CDim>
std::vector<
    LinearQuadraticStage<Scalar, Transcription::XDim, Transcription::UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<Scalar, Transcription, ConstraintConfig>&
        optimalControlProblem,
    const ContinuousTrajectory<Scalar, Transcription::XDim, Transcription::UDim,
                               Transcription::PredictLength>& nominalTrajectory,
    const StateInputConstraint<Scalar, Transcription::XDim, Transcription::UDim,
                               CDim>& constraint) {
  const auto targetTrajectory =
      toTargetTrajectories<Scalar, Transcription>(nominalTrajectory);
  std::array<MultiplierCollection<
                 Scalar, IntermediateStageConstraintLayout<ConstraintConfig>>,
             Transcription::PredictLength>
      intermediateMultipliers{};
  MultiplierCollection<Scalar, FinalStageConstraintLayout<ConstraintConfig>>
      finalMultipliers{};
  return getLinearQuadraticApproximation(
      optimalControlProblem, targetTrajectory, nominalTrajectory,
      intermediateMultipliers, finalMultipliers, constraint);
}

}  // namespace qp_solver
