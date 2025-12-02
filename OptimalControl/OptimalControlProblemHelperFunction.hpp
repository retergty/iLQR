#pragma once
#include "OptimalControlProblem.hpp"
#include "PrimalSolution.hpp"
#include "DualSolution.hpp"
#include "ProblemMetrics.hpp"
#include "Metrics.hpp"

/**
 * Initializes final MultiplierCollection for equality and inequality Lagrangians.
 *
 * @param [in] ocp : A const reference to the optimal control problem.
 * @param [in] time : Final time.
 * @param [out] multiplierCollection : The initialized final MultiplierCollection.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void initializeFinalMultiplierCollection(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                                         Scalar time, MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &multiplierCollection)
{
  ocp.finalEqualityLagrangian.initializeLagrangian(time, multiplierCollection.stateEq);
  ocp.finalInequalityLagrangian.initializeLagrangian(time, multiplierCollection.stateIneq);
}

/**
 * Initializes intermediate MultiplierCollection for equality and inequality Lagrangians.
 *
 * @param [in] ocp : A const reference to the optimal control problem.
 * @param [in] time : Intermediate time.
 * @param [out] multiplierCollection : The initialized intermediate MultiplierCollection.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void initializeIntermediateMultiplierCollection(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                                                Scalar time, MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multiplierCollection)
{
  ocp.stateEqualityLagrangian.initializeLagrangian(time, multiplierCollection.stateEq);
  ocp.stateInequalityLagrangian.initializeLagrangian(time, multiplierCollection.stateIneq);
  ocp.equalityLagrangian.initializeLagrangian(time, multiplierCollection.stateInputEq);
  ocp.inequalityLagrangian.initializeLagrangian(time, multiplierCollection.stateInputIneq);
}

/**
 * Initializes the dual solution based on the cached dual solution. It will use interpolation if cachedDualSolution has any component
 * in the same mode otherwise it will use the Lagrangian initialization method of ocp.
 *
 * @param [in] ocp : A const reference to the optimal control problem.
 * @param [in] primalSolution : The primal solution.
 * @param [in] cachedDualSolution : The cached dual solution which will be used for interpolation.
 * @param [out] dualSolution : The initialized dual solution.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void initializeDualSolution(
    const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
    const PrimalSolution<Scalar, XDimisions, UDimisions, PredictLength> &primalSolution,
    const DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength> &cachedDualSolution,
    DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength> &dualSolution)
{
  dualSolution.timeTrajectory = primalSolution.timeTrajectory_;

  if (!cachedDualSolution.empty())
  {
    // final
    dualSolution.final = cachedDualSolution.final;
  }
  else
  {
    initializeFinalMultiplierCollection(ocp, primalSolution.timeTrajectory_.back(), dualSolution.final);
  }

  if (!cachedDualSolution.empty())
  {
    // intermediates
    for (size_t i = 0; i < PredictLength; i++)
    {
      const Scalar &time = primalSolution.timeTrajectory_[i];
      MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multipliers = dualSolution.intermediates[i];
      multipliers = getIntermediateDualSolutionAtTime(cachedDualSolution, time);
    }
  }
  else
  {
    // intermediates
    for (size_t i = 0; i < PredictLength; i++)
    {
      const Scalar &time = primalSolution.timeTrajectory_[i];
      MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multipliers = dualSolution.intermediates[i];
      initializeIntermediateMultiplierCollection(ocp, time, multipliers);
    }
  }
}

/**
 * Updates in-place the dual solution based on its current solution and state-input values using the Lagrangian update method in ocp.
 * Moreover it also updates the penalties of ProblemMetrics based on the update of dual solution.
 *
 * @param [in] ocp : A const reference to the optimal control problem.
 * @param [in] primalSolution : The primal solution.
 * @param [in, out] problemMetrics : The problem metric. Its penalties will be updated based on the update of dualSolution.
 * @param [out] dualSolution : The updated dual solution.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void updateDualSolution(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                        const PrimalSolution<Scalar, XDimisions, UDimisions, PredictLength> &primalSolution,
                        ProblemMetrics<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &problemMetrics,
                        DualSolutionRef<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength> dualSolution)
{
  // final
  {
    const Scalar &time = primalSolution.timeTrajectory_.back();
    const Vector<Scalar, XDimisions> &state = primalSolution.stateTrajectory_.back();
    Metrics<Scalar, XDimisions, UDimisions, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &metrics = problemMetrics.final;
    MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &multipliers = dualSolution.final;
    updateFinalMultiplierCollection(ocp, time, state, metrics, multipliers);
  }

  // intermediates
  // static_assert(dualSolution.intermediates.size() == primalSolution.timeTrajectory_.size());
  // static_assert(problemMetrics.intermediates.size() == primalSolution.timeTrajectory_.size());

  for (size_t i = 0; i < PredictLength; i++)
  {
    const Scalar &time = primalSolution.timeTrajectory_[i];
    const Vector<Scalar, XDimisions> &state = primalSolution.stateTrajectory_[i];
    const Vector<Scalar, UDimisions> &input = primalSolution.inputTrajectory_[i];
    Metrics<Scalar, XDimisions, UDimisions, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &metrics = problemMetrics.intermediates[i];
    MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multipliers = dualSolution.intermediates[i];

    updateIntermediateMultiplierCollection(ocp, time, state, input, metrics, multipliers);
  }
}

/**
 * Updates in-place final MultiplierCollection for equality and inequality Lagrangians.
 * Moreover it also updates the penalties of Metrics based on the update of multipliers.
 *
 * @param [in] ocp : A const reference to the optimal control problem.
 * @param [in] time : Final time.
 * @param [in] state : Final state.
 * @param [in, out] metrics: The final Metrics. Its penalties will be updated based on the update of multiplierCollection.
 * @param [out] multipliers : The updated final MultiplierCollection.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void updateFinalMultiplierCollection(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                                     Scalar time, const Vector<Scalar, XDimisions> &state,
                                     Metrics<Scalar, XDimisions, UDimisions, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &metrics,
                                     MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0> &multipliers)
{
  ocp.finalEqualityLagrangian.updateLagrangian(time, state, metrics.stateEqLagrangian, multipliers.stateEq);
  ocp.finalInequalityLagrangian.updateLagrangian(time, state, metrics.stateIneqLagrangian, multipliers.stateIneq);
}

/**
 * Updates in-place intermediate MultiplierCollection for equality and inequality Lagrangians.
 * Moreover it also updates the penalties of Metrics based on the update of multipliers.
 *
 * @param [in] ocp : A const reference to the optimal control problem.
 * @param [in] time : Intermediate time.
 * @param [in] state : Intermediate state.
 * @param [in] input : Intermediate input.
 * @param [in, out] metrics: The intermediate Metrics. Its penalties will be updated based on the update of multiplierCollection.
 * @param [out] multipliers : The updated Intermediate MultiplierCollection.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
void updateIntermediateMultiplierCollection(const OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &ocp,
                                            Scalar time, const Vector<Scalar, XDimisions> &state, const Vector<Scalar, UDimisions> &input,
                                            Metrics<Scalar, XDimisions, UDimisions, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &metrics,
                                            MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> &multipliers)
{
  ocp.stateEqualityLagrangian.updateLagrangian(time, state, metrics.stateEqLagrangian, multipliers.stateEq);
  ocp.stateInequalityLagrangian.updateLagrangian(time, state, metrics.stateIneqLagrangian, multipliers.stateIneq);
  ocp.equalityLagrangian.updateLagrangian(time, state, input, metrics.stateInputEqLagrangian, multipliers.stateInputEq);
  ocp.inequalityLagrangian.updateLagrangian(time, state, input, metrics.stateInputIneqLagrangian, multipliers.stateInputIneq);
}
