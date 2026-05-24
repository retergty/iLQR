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
#include "SensitivityIntegrator.hpp"
#include "StateInputConstraint.hpp"
#include "iLQRDescriptor.hpp"

namespace qp_solver {

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpTranscriptionConfig_t =
    TranscriptionConfig<Dimensions<XDim, UDim>, Horizon<PredictLength>>;

using QpConstraintConfig_t = ConstraintConfig<>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpDescriptor_t =
    iLQRDescriptor<Scalar,
                   QpTranscriptionConfig_t<Scalar, XDim, UDim, PredictLength>,
                   QpConstraintConfig_t>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength>
using QpOptimalControlProblem_t = typename LinearQuadraticApproximator<
    QpDescriptor_t<Scalar, XDim, UDim, PredictLength>>::OptimalControlProblem_t;

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

template <typename Scalar, int XDim, int UDim>
typename LinearQuadraticStage<Scalar, XDim, UDim>::DynamicsApproximation_t
approximateDynamics(SystemDynamicsBase<Scalar, XDim, UDim>& system,
                    TrajectoryRef<Scalar, XDim, UDim> start, Scalar dt) {
  EK2DynamicsDiscretizer<Scalar, XDim, UDim> discretizer;
  auto discreteDynamics =
      discretizer.sensitivityDiscretize(system, start.t, start.x, start.u, dt);
  // Use deviation dynamics without an affine defect term.
  discreteDynamics.f.setZero();
  return discreteDynamics;
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
  using Approximator_t = LinearQuadraticApproximator<Descriptor_t>;
  using LqStage_t = LinearQuadraticStage<Scalar, XDim, UDim>;

  const auto modelData = Approximator_t::approximateIntermediateLQ(
      optimalControlProblem, targetTrajectory, start.t, start.x, start.u,
      multipliers);

  LqStage_t lqStage;
  const Scalar dt = end.t - start.t;

  lqStage.cost = modelData.cost;
  lqStage.cost *= dt;

  // Discrete LQ deviation dynamics: dx[k+1] = A dx[k] + B du[k].
  // Nominal trajectory defects are not included in this transcription.
  lqStage.dynamics =
      approximateDynamics(*optimalControlProblem.dynamicsPtr, start, dt);

  // In this project constraint penalties are folded into the approximated cost
  // through AugmentedLagrangian.
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
  using Approximator_t = LinearQuadraticApproximator<Descriptor_t>;

  const auto& t = nominalTrajectory.timeTrajectory;
  const auto& x = nominalTrajectory.stateTrajectory;
  const auto& u = nominalTrajectory.inputTrajectory;
  constexpr size_t N = PredictLength;

  // LinearQuadraticProblem with N+1 elements. Terminal stage lqp[N].dynamics is
  // ignored.
  std::vector<LinearQuadraticStage<Scalar, XDim, UDim>> lqp;
  lqp.reserve(N + 1);
  for (size_t k = 0; k < N; ++k) {  // Intermediate stages
    lqp.emplace_back(approximateStage(optimalControlProblem, targetTrajectory,
                                      {t[k], x[k], u[k]}, {t[k + 1], x[k + 1]},
                                      intermediateMultipliers[k]));
  }

  auto modelData = Approximator_t::approximateFinalLQ(
      optimalControlProblem, targetTrajectory, t[N], x[N], finalMultipliers);
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
  using Approximator_t = LinearQuadraticApproximator<Descriptor_t>;

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

  auto modelData = Approximator_t::approximateFinalLQ(
      optimalControlProblem, targetTrajectory, t[N], x[N], finalMultipliers);
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
