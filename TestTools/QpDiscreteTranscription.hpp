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
#include <vector>

#include "LinearQuadraticApproximator.hpp"
#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "QpSolverTypes.hpp"
#include "QpTrajectories.hpp"
#include "SensitivityIntegrator.hpp"
#include "iLQRDescriptor.hpp"

namespace qp_solver {

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using QpTranscriptionConfig_t =
    TranscriptionConfig<Dimensions<XDim, UDim>, Horizon<PredictLength>>;

template <int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using StateEqQpGroupLayout_t = std::conditional_t<
    StateEqLagrangianConstrainNumbers == 0, ConstraintGroupLayout<>,
    ConstraintGroupLayout<ConstraintTerm<StateEqLagrangianConstrainNumbers>>>;

template <int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using StateIneqQpGroupLayout_t = std::conditional_t<
    StateIneqLagrangianConstrainNumbers == 0, ConstraintGroupLayout<>,
    ConstraintGroupLayout<ConstraintTerm<StateIneqLagrangianConstrainNumbers>>>;

template <int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using StateInputEqQpGroupLayout_t = std::conditional_t<
    StateInputEqLagrangianConstrainNumbers == 0, ConstraintGroupLayout<>,
    ConstraintGroupLayout<
        ConstraintTerm<StateInputEqLagrangianConstrainNumbers>>>;

template <int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using StateInputIneqQpGroupLayout_t = std::conditional_t<
    StateInputIneqLagrangianConstrainNumbers == 0, ConstraintGroupLayout<>,
    ConstraintGroupLayout<
        ConstraintTerm<StateInputIneqLagrangianConstrainNumbers>>>;

template <int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using FinalStateEqQpGroupLayout_t = std::conditional_t<
    FinalStateEqLagrangianConstrainNumbers == 0, ConstraintGroupLayout<>,
    ConstraintGroupLayout<
        ConstraintTerm<FinalStateEqLagrangianConstrainNumbers>>>;

template <int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using FinalStateIneqQpGroupLayout_t = std::conditional_t<
    FinalStateIneqFinalLagrangianConstrainNumbers == 0, ConstraintGroupLayout<>,
    ConstraintGroupLayout<
        ConstraintTerm<FinalStateIneqFinalLagrangianConstrainNumbers>>>;

template <int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using QpConstraintConfig_t = ConstraintConfig<
    StateConstraintConfig<ConstraintLayout<
        StateEqQpGroupLayout_t<StateEqLagrangianConstrainNumbers,
                               StateIneqLagrangianConstrainNumbers,
                               StateInputEqLagrangianConstrainNumbers,
                               StateInputIneqLagrangianConstrainNumbers,
                               FinalStateEqLagrangianConstrainNumbers,
                               FinalStateIneqFinalLagrangianConstrainNumbers>,
        StateIneqQpGroupLayout_t<
            StateEqLagrangianConstrainNumbers,
            StateIneqLagrangianConstrainNumbers,
            StateInputEqLagrangianConstrainNumbers,
            StateInputIneqLagrangianConstrainNumbers,
            FinalStateEqLagrangianConstrainNumbers,
            FinalStateIneqFinalLagrangianConstrainNumbers>>>,
    StateInputConstraintConfig<
        ConstraintLayout<StateInputEqQpGroupLayout_t<
                             StateEqLagrangianConstrainNumbers,
                             StateIneqLagrangianConstrainNumbers,
                             StateInputEqLagrangianConstrainNumbers,
                             StateInputIneqLagrangianConstrainNumbers,
                             FinalStateEqLagrangianConstrainNumbers,
                             FinalStateIneqFinalLagrangianConstrainNumbers>,
                         StateInputIneqQpGroupLayout_t<
                             StateEqLagrangianConstrainNumbers,
                             StateIneqLagrangianConstrainNumbers,
                             StateInputEqLagrangianConstrainNumbers,
                             StateInputIneqLagrangianConstrainNumbers,
                             FinalStateEqLagrangianConstrainNumbers,
                             FinalStateIneqFinalLagrangianConstrainNumbers>>>,
    FinalStateConstraintConfig<
        ConstraintLayout<FinalStateEqQpGroupLayout_t<
                             StateEqLagrangianConstrainNumbers,
                             StateIneqLagrangianConstrainNumbers,
                             StateInputEqLagrangianConstrainNumbers,
                             StateInputIneqLagrangianConstrainNumbers,
                             FinalStateEqLagrangianConstrainNumbers,
                             FinalStateIneqFinalLagrangianConstrainNumbers>,
                         FinalStateIneqQpGroupLayout_t<
                             StateEqLagrangianConstrainNumbers,
                             StateIneqLagrangianConstrainNumbers,
                             StateInputEqLagrangianConstrainNumbers,
                             StateInputIneqLagrangianConstrainNumbers,
                             FinalStateEqLagrangianConstrainNumbers,
                             FinalStateIneqFinalLagrangianConstrainNumbers>>>>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using QpDescriptor_t = iLQRDescriptor<
    Scalar,
    QpTranscriptionConfig_t<Scalar, XDim, UDim, PredictLength,
                            StateEqLagrangianConstrainNumbers,
                            StateIneqLagrangianConstrainNumbers,
                            StateInputEqLagrangianConstrainNumbers,
                            StateInputIneqLagrangianConstrainNumbers,
                            FinalStateEqLagrangianConstrainNumbers,
                            FinalStateIneqFinalLagrangianConstrainNumbers>,
    QpConstraintConfig_t<StateEqLagrangianConstrainNumbers,
                         StateIneqLagrangianConstrainNumbers,
                         StateInputEqLagrangianConstrainNumbers,
                         StateInputIneqLagrangianConstrainNumbers,
                         FinalStateEqLagrangianConstrainNumbers,
                         FinalStateIneqFinalLagrangianConstrainNumbers>>;

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using QpOptimalControlProblem_t =
    typename LinearQuadraticApproximator<QpDescriptor_t<
        Scalar, XDim, UDim, PredictLength, StateEqLagrangianConstrainNumbers,
        StateIneqLagrangianConstrainNumbers,
        StateInputEqLagrangianConstrainNumbers,
        StateInputIneqLagrangianConstrainNumbers,
        FinalStateEqLagrangianConstrainNumbers,
        FinalStateIneqFinalLagrangianConstrainNumbers>>::
        OptimalControlProblem_t;

template <typename Scalar, int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers>
using QpIntermediateMultiplierCollection_t = MultiplierCollection<
    Scalar,
    IntermediateStageConstraintLayout<QpConstraintConfig_t<
        StateEqLagrangianConstrainNumbers, StateIneqLagrangianConstrainNumbers,
        StateInputEqLagrangianConstrainNumbers,
        StateInputIneqLagrangianConstrainNumbers, 0, 0>>>;

template <typename Scalar, int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
using QpFinalMultiplierCollection_t =
    MultiplierCollection<Scalar,
                         FinalStageConstraintLayout<QpConstraintConfig_t<
                             0, 0, 0, 0, FinalStateEqLagrangianConstrainNumbers,
                             FinalStateIneqFinalLagrangianConstrainNumbers>>>;

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
      optimalControlProblem, start.t, start.x, start.u, multipliers);

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
    lqp.emplace_back(approximateStage(optimalControlProblem, {t[k], x[k], u[k]},
                                      {t[k + 1], x[k + 1]},
                                      intermediateMultipliers[k]));
  }

  auto modelData = Approximator_t::approximateFinalLQ(
      optimalControlProblem, t[N], x[N], finalMultipliers);
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
      optimalControlProblem, nominalTrajectory, intermediateMultipliers,
      finalMultipliers);
}

}  // namespace qp_solver
