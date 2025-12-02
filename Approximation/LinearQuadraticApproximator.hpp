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

#include "Types.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "OptimalControlProblem.hpp"
#include "Multiplier.hpp"
#include "ModelData.hpp"
#include "Metrics.hpp"

template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqLagrangianConstrains, int StateIneqLagrangianConstrains, int StateInputEqLagrangianConstrains, int StateInputIneqLagrangianConstrains,
          int FinalStateEqLagrangianConstrains, int FinalStateIneqFinalLagrangianConstrains>
struct LinearQuadraticApproximator
{
  using OptimalControlProblem_t = OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqLagrangianConstrains, StateInputEqLagrangianConstrains, StateIneqLagrangianConstrains, StateInputIneqLagrangianConstrains, FinalStateEqLagrangianConstrains, FinalStateIneqFinalLagrangianConstrains>;
  using StateVector_t = Vector<Scalar, XDimisions>;
  using InputVector_t = Vector<Scalar, UDimisions>;
  using ModelData_t = ModelData<Scalar, XDimisions, UDimisions>;
  using IntermediateMultiplierCollection_t = MultiplierCollection<Scalar, StateEqLagrangianConstrains, StateIneqLagrangianConstrains, StateInputEqLagrangianConstrains, StateInputIneqLagrangianConstrains>;
  using FinalMultiplierCollection_t = MultiplierCollection<Scalar, FinalStateEqLagrangianConstrains, FinalStateIneqFinalLagrangianConstrains, 0, 0>;
  using IntermediateMetrics_t = Metrics<Scalar, XDimisions, UDimisions, StateEqLagrangianConstrains, StateIneqLagrangianConstrains, StateInputEqLagrangianConstrains, StateInputIneqLagrangianConstrains>;
  using FinalMetrics_t = Metrics<Scalar, XDimisions, UDimisions, FinalStateEqLagrangianConstrains, FinalStateIneqFinalLagrangianConstrains, 0, 0>;
  using TimeTrajectory_t = std::array<Scalar, PredictLength + 1>;
  using StateTrajectory_t = std::array<Vector<Scalar, XDimisions>, PredictLength + 1>;
  using InputTrajectory_t = std::array<Vector<Scalar, UDimisions>, PredictLength + 1>;

  /**
   * Calculates an LQ approximate of the constrained optimal control problem at a given time, state, and input.
   *
   * @param [in] problem: The optimal control problem
   * @param [in] time: The current time.
   * @param [in] state: The current state.
   * @param [in] input: The current input.
   * @param [in] multipliers: The current multipliers associated to the equality and inequality Lagrangians.
   * @param [out] modelData: The output data model.
   */
  static void approximateIntermediateLQ(const OptimalControlProblem_t &problem, const Scalar time, const StateVector_t &state, const InputVector_t &input,
                                        const IntermediateMultiplierCollection_t &multipliers, ModelData_t &modelData)
  {
    // Dynamics
    modelData.dynamics = problem.dynamicsPtr->linearApproximation(time, state, input);

    // Cost
    modelData.cost = approximateCost(problem, time, state, input);

    // Lagrangians
    if constexpr (StateEqLagrangianConstrains != 0)
    {
      ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> approx = problem.stateEqualityLagrangian.getQuadraticApproximation(time, state, multipliers.stateEq);
      modelData.cost.f += approx.f;
      modelData.cost.dfdx += approx.dfdx;
      modelData.cost.dfdxx += approx.dfdxx;
    }
    if constexpr (StateIneqLagrangianConstrains != 0)
    {
      ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> approx = problem.stateInequalityLagrangian.getQuadraticApproximation(time, state, multipliers.stateIneq);
      modelData.cost.f += approx.f;
      modelData.cost.dfdx += approx.dfdx;
      modelData.cost.dfdxx += approx.dfdxx;
    }
    if constexpr (StateInputEqLagrangianConstrains != 0)
    {
      modelData.cost +=
          problem.equalityLagrangian.getQuadraticApproximation(time, state, input, multipliers.stateInputEq);
    }
    if constexpr (StateInputIneqLagrangianConstrains != 0)
    {
      modelData.cost +=
          problem.inequalityLagrangian.getQuadraticApproximation(time, state, input, multipliers.stateInputIneq);
    }
  }

  /**
   * Calculates an LQ approximate of the constrained optimal control problem at a given time, state, and input.
   *
   * @param [in] problem: The optimal control problem
   * @param [in] time: The current time.
   * @param [in] state: The current state.
   * @param [in] input: The current input.
   * @param [in] multipliers: The current multipliers associated to the equality and inequality Lagrangians.
   * @return The output data model.
   */
  static inline ModelData_t approximateIntermediateLQ(const OptimalControlProblem_t &problem, const Scalar time, const StateVector_t &state, const InputVector_t &input,
                                                      const IntermediateMultiplierCollection_t &multipliers)
  {
    ModelData<Scalar, XDimisions, UDimisions> md;
    approximateIntermediateLQ(problem, time, state, input, multipliers, md);
    return md;
  }

  /**
   * Calculates an LQ approximate of the constrained optimal control problem at final time.
   *
   * @param [in] problem: The optimal control problem
   * @param [in] time: The current time.
   * @param [in] state: The current state.
   * @param [in] multipliers: The current multipliers associated to the equality and inequality Lagrangians.
   * @param [out] modelData: The output data model.
   */
  static void approximateFinalLQ(const OptimalControlProblem_t &problem,
                                 const Scalar time, const StateVector_t &state,
                                 const FinalMultiplierCollection_t &multipliers, ModelData_t &modelData)
  {
    modelData.time = time;

    VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions> finalDynamics;
    finalDynamics.setZero();

    // Dynamics
    modelData.dynamics = finalDynamics;

    // Final cost
    modelData.cost = approximateFinalCost(problem, time, state);

    // Lagrangians
    if constexpr (FinalStateEqLagrangianConstrains != 0)
    {
      auto approx = problem.finalEqualityLagrangian.getQuadraticApproximation(time, state, multipliers.stateEq);
      modelData.cost.f += approx.f;
      modelData.cost.dfdx += approx.dfdx;
      modelData.cost.dfdxx += approx.dfdxx;
    }
    if constexpr (FinalStateIneqFinalLagrangianConstrains != 0)
    {
      auto approx = problem.finalInequalityLagrangian.getQuadraticApproximation(time, state, multipliers.stateIneq);
      modelData.cost.f += approx.f;
      modelData.cost.dfdx += approx.dfdx;
      modelData.cost.dfdxx += approx.dfdxx;
    }
  }

  /**
   * Calculates an LQ approximate of the constrained optimal control problem at final time.
   *
   * @param [in] problem: The optimal control problem
   * @param [in] time: The current time.
   * @param [in] state: The current state.
   * @param [in] multipliers: The current multipliers associated to the equality and inequality Lagrangians.
   * @return The output data model.
   */
  static inline ModelData_t approximateFinalLQ(const OptimalControlProblem_t &problem,
                                               const Scalar time, const StateVector_t &state,
                                               const FinalMultiplierCollection_t &multipliers)
  {
    ModelData_t md;
    approximateFinalLQ(problem, time, state, multipliers, md);
    return md;
  }

  /**
   * Compute the total intermediate cost (i.e. cost + softConstraints). It is assumed that the precomputation request is already made.
   */
  static Scalar computeCost(const OptimalControlProblem_t &problem, const Scalar time, const StateVector_t &state, const InputVector_t &input)
  {
    const TimeTrajectory_t &targetTimeTrajectories = problem.timeTrajectory;
    const StateTrajectory_t &targetStateTrajectories = problem.stateTrajectory;
    const InputTrajectory_t &targetInputTrajectories = problem.inputTrajectory;

    // Compute and sum all costs
    Scalar cost = problem.cost.getValue(time, state, input, targetTimeTrajectories, targetStateTrajectories, targetInputTrajectories);
    cost += problem.stateCost.getValue(time, state, targetTimeTrajectories, targetStateTrajectories);

    return cost;
  }

  /**
   * Compute the quadratic approximation of the total intermediate cost (i.e. cost + softConstraints). It is assumed that the precomputation
   * request is already made.
   */
  static ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>
  approximateCost(const OptimalControlProblem_t &problem,
                  const Scalar time, const StateVector_t &state, const InputVector_t &input)
  {
    const TimeTrajectory_t &targetTimeTrajectories = problem.timeTrajectory;
    const StateTrajectory_t &targetStateTrajectories = problem.stateTrajectory;
    const InputTrajectory_t &targetInputTrajectories = problem.inputTrajectory;

    // get the state-input cost approximations
    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> cost = problem.cost.getQuadraticApproximation(time, state, input, targetTimeTrajectories, targetStateTrajectories, targetInputTrajectories);

    // get the state only cost approximations
    cost += problem.stateCost.getQuadraticApproximation(time, state, targetTimeTrajectories, targetStateTrajectories);

    return cost;
  }

  /**
   * Compute the total final cost (i.e. cost + softConstraints). It is assumed that the precomputation request is already made.
   */
  static Scalar computeFinalCost(const OptimalControlProblem_t &problem,
                                 const Scalar time, const StateVector_t &state)
  {
    const TimeTrajectory_t &targetTimeTrajectories = problem.timeTrajectory;
    const StateTrajectory_t &targetStateTrajectories = problem.stateTrajectory;

    Scalar cost = problem.finalCost.getValue(time, state, targetTimeTrajectories, targetStateTrajectories);

    return cost;
  }

  /**
   * Compute the quadratic approximation of the total final cost (i.e. cost + softConstraints). It is assumed that the precomputation
   * request is already made.
   */
  static ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>
  approximateFinalCost(const OptimalControlProblem_t &problem,
                       const Scalar time, const StateVector_t &state)
  {
    const TimeTrajectory_t &targetTimeTrajectories = problem.timeTrajectory;
    const StateTrajectory_t &targetStateTrajectories = problem.stateTrajectory;

    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> cost = problem.finalCost.getQuadraticApproximation(time, state, targetTimeTrajectories, targetStateTrajectories);

    return cost;
  }

  /**
   * Compute the intermediate-time Metrics (i.e. cost, softConstraints, and constraints).
   *
   * @note It is assumed that the precomputation request is already made.
   * problem.preComputationPtr->request(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x, u)
   *
   * @param [in] problem: The optimal control probelm
   * @param [in] time: The current time.
   * @param [in] state: The current state.
   * @param [in] input: The current input.
   * @param [in] dynamicsViolation: The violation of dynamics. It depends on the transcription method.
   * @return The output Metrics.
   */
  static IntermediateMetrics_t
  computeIntermediateMetrics(const OptimalControlProblem_t &problem,
                             const Scalar time, const StateVector_t &state, const InputVector_t &input)
  {
    IntermediateMetrics_t metrics;

    // Cost
    metrics.cost = computeCost(problem, time, state, input);

    return metrics;
  }

  /**
   * Compute the intermediate-time Metrics (i.e. cost, softConstraints, and constraints).
   *
   * @note It is assumed that the precomputation request is already made.
   * problem.preComputationPtr->request(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x, u)
   *
   * @param [in] problem: The optimal control probelm
   * @param [in] time: The current time.
   * @param [in] state: The current state.
   * @param [in] input: The current input.
   * @param [in] multipliers: The current multipliers associated to the equality and inequality Lagrangians.
   * @param [in] dynamicsViolation: The violation of dynamics. It depends on the transcription method.
   * @return The output Metrics.
   */
  static IntermediateMetrics_t
  computeIntermediateMetrics(const OptimalControlProblem_t &problem,
                             const Scalar time, const StateVector_t &state, const InputVector_t &input,
                             const IntermediateMultiplierCollection_t &multipliers)
  {
    // cost, dynamics violation, equlaity constraints, inequlaity constraints
    IntermediateMetrics_t metrics = computeIntermediateMetrics(problem, time, state, input);

    // Equality Lagrangians
    if constexpr (StateEqLagrangianConstrains != 0)
    {
      metrics.stateEqLagrangian = problem.stateEqualityLagrangian.getValue(time, state, multipliers.stateEq);
    }

    if constexpr (StateInputEqLagrangianConstrains != 0)
    {
      metrics.stateInputEqLagrangian = problem.equalityLagrangian.getValue(time, state, input, multipliers.stateInputEq);
    }

    // Inequality Lagrangians
    if constexpr (StateIneqLagrangianConstrains != 0)
    {
      metrics.stateIneqLagrangian = problem.stateInequalityLagrangian.getValue(time, state, multipliers.stateIneq);
    }

    if constexpr (StateInputIneqLagrangianConstrains != 0)
    {
      metrics.stateInputIneqLagrangian = problem.inequalityLagrangian.getValue(time, state, input, multipliers.stateInputIneq);
    }

    return metrics;
  }

  /**
   * Compute the final-time Metrics (i.e. cost, softConstraints, and constraints).
   *
   * @note It is assumed that the precomputation request is already made.
   * problem.preComputationPtr->requestFinal(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x)
   *
   * @param [in] problem: The optimal control probelm
   * @param [in] time: The current time.
   * @param [in] state: The current state.
   * @return The output Metrics.
   */
  static FinalMetrics_t
  computeFinalMetrics(const OptimalControlProblem_t &problem,
                      const Scalar time, const StateVector_t &state)
  {

    FinalMetrics_t metrics;

    // Cost
    metrics.cost = computeFinalCost(problem, time, state);

    return metrics;
  }

  /**
   * Compute the final-time Metrics (i.e. cost, softConstraints, and constraints).
   *
   * @note It is assumed that the precomputation request is already made.
   * problem.preComputationPtr->requestFinal(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x)
   *
   * @param [in] problem: The optimal control probelm
   * @param [in] time: The current time.
   * @param [in] state: The current state.
   * @param [in] multipliers: The current multipliers associated to the equality and inequality Lagrangians.
   * @return The output Metrics.
   */
  static FinalMetrics_t
  computeFinalMetrics(const OptimalControlProblem_t &problem,
                      const Scalar time, const StateVector_t &state, const FinalMultiplierCollection_t &multipliers)
  {

    // cost, equlaity constraints, inequlaity constraints
    FinalMetrics_t metrics = computeFinalMetrics(problem, time, state);

    if constexpr (FinalStateEqLagrangianConstrains != 0)
    {
      // Equality Lagrangians
      metrics.stateEqLagrangian = problem.finalEqualityLagrangian.getValue(time, state, multipliers.stateEq);
    }

    // Inequality Lagrangians
    if constexpr (FinalStateIneqFinalLagrangianConstrains != 0)
    {
      metrics.stateIneqLagrangian = problem.finalInequalityLagrangian.getValue(time, state, multipliers.stateIneq);
    }

    return metrics;
  }
};
