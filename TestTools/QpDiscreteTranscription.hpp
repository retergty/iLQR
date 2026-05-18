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
#include "ModelData.hpp"
#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "QpSolverTypes.hpp"
#include "QpTrajectories.hpp"

namespace qp_solver {

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
approximateDynamics(const ModelData<Scalar, XDim, UDim>& modelData,
                    TrajectoryRef<Scalar, XDim, UDim> start, Scalar dt) {
  // Forward Euler discretization
  // x[k+1] = x[k] + dt * dxdt[k]
  // x[k+1] = (x0[k] + dx[k]) + dt * dxdt[k]
  // x[k+1] = (x0[k] + dx[k]) + dt * (A_c dx[k] + B_c du[k] + b_c)
  // x[k+1] = (I + A_c * dt) dx[k] + (B_c * dt) du[k] + (b_c * dt + x0[k])
  const auto& continuousDynamics = modelData.dynamics;
  typename LinearQuadraticStage<Scalar, XDim, UDim>::DynamicsApproximation_t
      discreteDynamics;
  discreteDynamics.dfdx = continuousDynamics.dfdx * dt;
  discreteDynamics.dfdx.diagonal().array() += 1.0;
  discreteDynamics.dfdu = continuousDynamics.dfdu * dt;
  discreteDynamics.f = continuousDynamics.f * dt + start.x;
  return discreteDynamics;
}

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
LinearQuadraticStage<Scalar, XDim, UDim> approximateStage(
    const OptimalControlProblem<
        Scalar, XDim, UDim, PredictLength, StateEqLagrangianConstrainNumbers,
        StateInputEqLagrangianConstrainNumbers,
        StateIneqLagrangianConstrainNumbers,
        StateInputIneqLagrangianConstrainNumbers,
        FinalStateEqLagrangianConstrainNumbers,
        FinalStateIneqFinalLagrangianConstrainNumbers>& optimalControlProblem,
    TrajectoryRef<Scalar, XDim, UDim> start,
    StateTrajectoryRef<Scalar, XDim> end,
    const MultiplierCollection<Scalar, StateEqLagrangianConstrainNumbers,
                               StateIneqLagrangianConstrainNumbers,
                               StateInputEqLagrangianConstrainNumbers,
                               StateInputIneqLagrangianConstrainNumbers>&
        multipliers) {
  using Approximator_t = LinearQuadraticApproximator<
      Scalar, XDim, UDim, PredictLength, StateEqLagrangianConstrainNumbers,
      StateIneqLagrangianConstrainNumbers,
      StateInputEqLagrangianConstrainNumbers,
      StateInputIneqLagrangianConstrainNumbers,
      FinalStateEqLagrangianConstrainNumbers,
      FinalStateIneqFinalLagrangianConstrainNumbers>;
  using LqStage_t = LinearQuadraticStage<Scalar, XDim, UDim>;

  const auto modelData = Approximator_t::approximateIntermediateLQ(
      optimalControlProblem, start.t, start.x, start.u, multipliers);

  LqStage_t lqStage;
  const Scalar dt = end.t - start.t;

  lqStage.cost = modelData.cost;
  lqStage.cost *= dt;

  // Linearized Dynamics after discretization: x0[k+1] + dx[k+1] = A dx[k] + B
  // du[k] + F(x0[k], u0[k])
  lqStage.dynamics = approximateDynamics(modelData, start, dt);
  // Adapt the offset to account for discretization and the nominal trajectory:
  // dx[k+1] = A dx[k] + B du[k] + F(x0[k], u0[k]) - x0[k+1]
  lqStage.dynamics.f -= end.x;

  // In this project constraint penalties are folded into the approximated cost
  // through AugmentedLagrangian.
  lqStage.constraints = makeZeroConstraints<Scalar, XDim, UDim>();

  return lqStage;
}

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<
        Scalar, XDim, UDim, PredictLength, StateEqLagrangianConstrainNumbers,
        StateInputEqLagrangianConstrainNumbers,
        StateIneqLagrangianConstrainNumbers,
        StateInputIneqLagrangianConstrainNumbers,
        FinalStateEqLagrangianConstrainNumbers,
        FinalStateIneqFinalLagrangianConstrainNumbers>& optimalControlProblem,
    const ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>&
        nominalTrajectory,
    const std::array<
        MultiplierCollection<Scalar, StateEqLagrangianConstrainNumbers,
                             StateIneqLagrangianConstrainNumbers,
                             StateInputEqLagrangianConstrainNumbers,
                             StateInputIneqLagrangianConstrainNumbers>,
        PredictLength>& intermediateMultipliers,
    const MultiplierCollection<Scalar, FinalStateEqLagrangianConstrainNumbers,
                               FinalStateIneqFinalLagrangianConstrainNumbers, 0,
                               0>& finalMultipliers) {
  if (optimalControlProblem.dynamicsPtr == nullptr) {
    throw std::runtime_error(
        "[qp_solver::getLinearQuadraticApproximation] "
        "dynamicsPtr should not be null.");
  }

  using Approximator_t = LinearQuadraticApproximator<
      Scalar, XDim, UDim, PredictLength, StateEqLagrangianConstrainNumbers,
      StateIneqLagrangianConstrainNumbers,
      StateInputEqLagrangianConstrainNumbers,
      StateInputIneqLagrangianConstrainNumbers,
      FinalStateEqLagrangianConstrainNumbers,
      FinalStateIneqFinalLagrangianConstrainNumbers>;

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

template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>
getLinearQuadraticApproximation(
    const OptimalControlProblem<
        Scalar, XDim, UDim, PredictLength, StateEqLagrangianConstrainNumbers,
        StateInputEqLagrangianConstrainNumbers,
        StateIneqLagrangianConstrainNumbers,
        StateInputIneqLagrangianConstrainNumbers,
        FinalStateEqLagrangianConstrainNumbers,
        FinalStateIneqFinalLagrangianConstrainNumbers>& optimalControlProblem,
    const ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>&
        nominalTrajectory) {
  std::array<MultiplierCollection<Scalar, StateEqLagrangianConstrainNumbers,
                                  StateIneqLagrangianConstrainNumbers,
                                  StateInputEqLagrangianConstrainNumbers,
                                  StateInputIneqLagrangianConstrainNumbers>,
             PredictLength>
      intermediateMultipliers{};
  MultiplierCollection<Scalar, FinalStateEqLagrangianConstrainNumbers,
                       FinalStateIneqFinalLagrangianConstrainNumbers, 0, 0>
      finalMultipliers{};
  return getLinearQuadraticApproximation(
      optimalControlProblem, nominalTrajectory, intermediateMultipliers,
      finalMultipliers);
}

}  // namespace qp_solver
