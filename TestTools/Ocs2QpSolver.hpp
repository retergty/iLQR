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
#include <utility>

#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "QpDiscreteTranscription.hpp"
#include "QpSolver.hpp"
#include "QpTrajectories.hpp"

namespace qp_solver {

/**
 * Solves a constrained discrete-time linear quadratic control problem around a
 * provided linearization trajectory. The time horizon and discretization steps
 * are defined by the time trajectory of the provided linearization.
 *
 * @param optimalControlProblem: The optimal control problem definition.
 * @param nominalTrajectory : time, state and input trajectory to make the
 * linear quadratic approximation around
 * @param initialState : state at the start of the horizon.
 * @param intermediateMultipliers : multipliers for intermediate augmented
 * Lagrangians.
 * @param finalMultipliers : multipliers for final augmented Lagrangians.
 * @return time, state, and input solution.
 */
template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>
solveLinearQuadraticOptimalControlProblem(
    const OptimalControlProblem<
        Scalar, XDim, UDim, PredictLength, StateEqLagrangianConstrainNumbers,
        StateIneqLagrangianConstrainNumbers,
        StateInputEqLagrangianConstrainNumbers,
        StateInputIneqLagrangianConstrainNumbers,
        FinalStateEqLagrangianConstrainNumbers,
        FinalStateIneqFinalLagrangianConstrainNumbers>& optimalControlProblem,
    const ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>&
        nominalTrajectory,
    const Vector<Scalar, XDim>& initialState,
    const std::array<
        MultiplierCollection<Scalar, StateEqLagrangianConstrainNumbers,
                             StateIneqLagrangianConstrainNumbers,
                             StateInputEqLagrangianConstrainNumbers,
                             StateInputIneqLagrangianConstrainNumbers>,
        PredictLength>& intermediateMultipliers,
    const MultiplierCollection<Scalar, FinalStateEqLagrangianConstrainNumbers,
                               FinalStateIneqFinalLagrangianConstrainNumbers, 0,
                               0>& finalMultipliers) {
  // Approximate
  const auto lqApproximation = getLinearQuadraticApproximation(
      optimalControlProblem, nominalTrajectory, intermediateMultipliers,
      finalMultipliers);

  // Solve for an update step
  ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> deltaSolution;
  deltaSolution.timeTrajectory = nominalTrajectory.timeTrajectory;
  const Vector<Scalar, XDim> dx0 =
      initialState - nominalTrajectory.stateTrajectory.front();
  const auto deltaTrajectories =
      solveLinearQuadraticProblem(lqApproximation, dx0);
  const auto& stateDeltaTrajectory = deltaTrajectories.first;
  const auto& inputDeltaTrajectory = deltaTrajectories.second;

  for (size_t k = 0; k < PredictLength; ++k) {
    deltaSolution.stateTrajectory[k] = stateDeltaTrajectory[k];
    deltaSolution.inputTrajectory[k] = inputDeltaTrajectory[k];
  }
  deltaSolution.stateTrajectory[PredictLength] =
      stateDeltaTrajectory[PredictLength];

  // Take a full step: Add update to nominal trajectory
  return nominalTrajectory + deltaSolution;
}

/**
 * Overload for problems without an externally provided multiplier trajectory.
 */
template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqLagrangianConstrainNumbers,
          int StateIneqLagrangianConstrainNumbers,
          int StateInputEqLagrangianConstrainNumbers,
          int StateInputIneqLagrangianConstrainNumbers,
          int FinalStateEqLagrangianConstrainNumbers,
          int FinalStateIneqFinalLagrangianConstrainNumbers>
ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>
solveLinearQuadraticOptimalControlProblem(
    const OptimalControlProblem<
        Scalar, XDim, UDim, PredictLength, StateEqLagrangianConstrainNumbers,
        StateIneqLagrangianConstrainNumbers,
        StateInputEqLagrangianConstrainNumbers,
        StateInputIneqLagrangianConstrainNumbers,
        FinalStateEqLagrangianConstrainNumbers,
        FinalStateIneqFinalLagrangianConstrainNumbers>& optimalControlProblem,
    const ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>&
        nominalTrajectory,
    const Vector<Scalar, XDim>& initialState) {
  std::array<MultiplierCollection<Scalar, StateEqLagrangianConstrainNumbers,
                                  StateIneqLagrangianConstrainNumbers,
                                  StateInputEqLagrangianConstrainNumbers,
                                  StateInputIneqLagrangianConstrainNumbers>,
             PredictLength>
      intermediateMultipliers{};
  MultiplierCollection<Scalar, FinalStateEqLagrangianConstrainNumbers,
                       FinalStateIneqFinalLagrangianConstrainNumbers, 0, 0>
      finalMultipliers{};
  return solveLinearQuadraticOptimalControlProblem(
      optimalControlProblem, nominalTrajectory, initialState,
      intermediateMultipliers, finalMultipliers);
}

}  // namespace qp_solver
