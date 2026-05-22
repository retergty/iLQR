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

#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "LinearQuadraticApproximator.hpp"
#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "QpSolverTypes.hpp"
#include "QpTrajectories.hpp"
#include "Reference.hpp"
#include "SensitivityIntegrator.hpp"
#include "iLQRDescriptor.hpp"

namespace qp_solver {

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqConstraintDim,
          int StateIneqConstraintDim,
          int StateInputEqConstraintDim,
          int StateInputIneqConstraintDim,
          int FinalStateEqConstraintDim,
          int FinalStateIneqConstraintDim>
using QpTranscriptionConfig_t =
    TranscriptionConfig<Dimensions<XDim, UDim>, Horizon<PredictLength>>;

template <int ConstraintDim>
using QpConstraintGroupLayout_t =
    std::conditional_t<ConstraintDim == 0, ConstraintGroupLayout<>,
                       ConstraintGroupLayout<ConstraintTerm<ConstraintDim>>>;

template <int StateEqConstraintDim,
          int StateIneqConstraintDim,
          int StateInputEqConstraintDim,
          int StateInputIneqConstraintDim,
          int FinalStateEqConstraintDim,
          int FinalStateIneqConstraintDim>
using QpConstraintConfig_t = ConstraintConfig<
    StateConstraintConfig<ConstraintLayout<
        QpConstraintGroupLayout_t<StateEqConstraintDim>,
        QpConstraintGroupLayout_t<StateIneqConstraintDim>>>,
    StateInputConstraintConfig<
        ConstraintLayout<QpConstraintGroupLayout_t<StateInputEqConstraintDim>,
                         QpConstraintGroupLayout_t<StateInputIneqConstraintDim>>>,
    FinalStateConstraintConfig<
        ConstraintLayout<QpConstraintGroupLayout_t<FinalStateEqConstraintDim>,
                         QpConstraintGroupLayout_t<FinalStateIneqConstraintDim>>>>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqConstraintDim,
          int StateIneqConstraintDim,
          int StateInputEqConstraintDim,
          int StateInputIneqConstraintDim,
          int FinalStateEqConstraintDim,
          int FinalStateIneqConstraintDim>
using QpDescriptor_t = iLQRDescriptor<
    Scalar,
    QpTranscriptionConfig_t<Scalar, XDim, UDim, PredictLength,
                            StateEqConstraintDim,
                            StateIneqConstraintDim,
                            StateInputEqConstraintDim,
                            StateInputIneqConstraintDim,
                            FinalStateEqConstraintDim,
                            FinalStateIneqConstraintDim>,
    QpConstraintConfig_t<StateEqConstraintDim,
                         StateIneqConstraintDim,
                         StateInputEqConstraintDim,
                         StateInputIneqConstraintDim,
                         FinalStateEqConstraintDim,
                         FinalStateIneqConstraintDim>>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqConstraintDim,
          int StateIneqConstraintDim,
          int StateInputEqConstraintDim,
          int StateInputIneqConstraintDim,
          int FinalStateEqConstraintDim,
          int FinalStateIneqConstraintDim>
using QpOptimalControlProblem_t =
    typename LinearQuadraticApproximator<QpDescriptor_t<
        Scalar, XDim, UDim, PredictLength, StateEqConstraintDim,
        StateIneqConstraintDim,
        StateInputEqConstraintDim,
        StateInputIneqConstraintDim,
        FinalStateEqConstraintDim,
        FinalStateIneqConstraintDim>>::
        OptimalControlProblem_t;

template <typename Scalar, int StateEqConstraintDim,
          int StateIneqConstraintDim,
          int StateInputEqConstraintDim,
          int StateInputIneqConstraintDim>
using QpIntermediateMultiplierCollection_t = MultiplierCollection<
    Scalar,
    IntermediateStageConstraintLayout<QpConstraintConfig_t<
        StateEqConstraintDim, StateIneqConstraintDim,
        StateInputEqConstraintDim,
        StateInputIneqConstraintDim, 0, 0>>>;

template <typename Scalar, int FinalStateEqConstraintDim,
          int FinalStateIneqConstraintDim>
using QpFinalMultiplierCollection_t =
    MultiplierCollection<Scalar,
                         FinalStageConstraintLayout<QpConstraintConfig_t<
                             0, 0, 0, 0, FinalStateEqConstraintDim,
                             FinalStateIneqConstraintDim>>>;

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
                                      {t[k], x[k], u[k]},
                                      {t[k + 1], x[k + 1]},
                                      intermediateMultipliers[k]));
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

}  // namespace qp_solver
