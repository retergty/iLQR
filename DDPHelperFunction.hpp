#pragma once

#include "Types.hpp"
#include "PrimalSolution.hpp"
#include "DualSolution.hpp"
#include "RolloutBase.hpp"
#include "TrapezoidalIntegration.hpp"
#include "PerformanceIndex.hpp"
/**
 * Forward integrate the system dynamics with given controller. It uses the given control policies and initial state,
 * to integrate the system dynamics in time period [initTime, finalTime].
 *
 * @param [in] rollout: A reference to the rollout class.
 * @param [in] initTime: The initial time.
 * @param [in] initState: The initial state.
 * @param [in] finalTime: The final time.
 * @param [in, out] primalSolution: The resulting primal solution. Make sure that primalSolution::controllerPtr is set since
 *                                  the rollout is performed based on the controller stored in primalSolution. Moreover,
 *                                  except for StateTriggeredRollout, one should also set primalSolution::modeSchedule.
 *
 * @return average time step.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength>
Scalar rolloutTrajectory(RolloutBase<Scalar, XDimisions, UDimisions, PredictLength> &rollout,
                         Scalar initTime, const Vector<Scalar, XDimisions> &initState, Scalar finalTime,
                         PrimalSolution<Scalar, XDimisions, UDimisions, PredictLength> &primalSolution)
{
  // rollout with controller
  // const Vector<Scalar, XDimisions> xCurrent = rollout.run(initTime, initState, finalTime, primalSolution.controllerPtr_,
  //                                                         primalSolution.timeTrajectory_, primalSolution.stateTrajectory_, primalSolution.inputTrajectory_);

  // assert(!xCurrent.allFinite());
  rollout.run(initTime, initState, finalTime, primalSolution.controllerPtr_,
              primalSolution.timeTrajectory_, primalSolution.stateTrajectory_, primalSolution.inputTrajectory_);
  // average time step
  return (finalTime - initTime) / static_cast<Scalar>(PredictLength);
}

/**
 * Calculates the PerformanceIndex associated to the given ProblemMetrics.
 *
 * @param [in] timeTrajectory: Time stamp of the rollout.
 * @param [in] problemMetrics: The cost, soft constraints and constraints values of the rollout.
 *
 * @return The PerformanceIndex of the trajectory.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
PerformanceIndex<Scalar> computeRolloutPerformanceIndex(
    const std::array<Scalar, PredictLength> &timeTrajectory,
    const ProblemMetrics<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains> &problemMetrics)
{
  // Final
  PerformanceIndex<Scalar> performanceIndex = toPerformanceIndex(problemMetrics.final);
  // Intermediates
  std::array<PerformanceIndex<Scalar>, PredictLength> performanceIndexTrajectory;

  // std::for_each(problemMetrics.intermediates.cbegin(), problemMetrics.intermediates.cend(),
  //   [&performanceIndexTrajectory](const Metrics& m) { performanceIndexTrajectory.push_back(toPerformanceIndex(m)); });

  for (size_t i = 0; i < PredictLength; ++i)
  {
    performanceIndexTrajectory[i] = toPerformanceIndex(problemMetrics.intermediates[i]);
  }

  // Intermediates
  return trapezoidalIntegration(timeTrajectory, performanceIndexTrajectory, performanceIndex);
}
/**
 * Outputs a controller with the same time stamp and gains as unoptimizedController. However, bias is incremented based on:
 * biasArray = unoptimizedController.biasArray + stepLength * unoptimizedController.deltaBiasArray
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength>
void incrementController(Scalar stepLength, const LinearController<Scalar, XDimisions, UDimisions, PredictLength> &unoptimizedController, LinearController<Scalar, XDimisions, UDimisions, PredictLength> &controller)
{
  controller.timeStamp_ = unoptimizedController.timeStamp_;
  controller.gainArray_ = unoptimizedController.gainArray_;
  for (size_t k = 0; k < PredictLength; k++)
  {
    controller.biasArray_[k] = unoptimizedController.biasArray_[k] + stepLength * unoptimizedController.deltaBiasArray_[k];
  }
}

/**
 * Computes the integral of the squared (IS) norm of the controller update.
 *
 * @param [in] controller: Input controller.
 * @return The integral of the squared (IS) norm of the controller update.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength>
Scalar computeControllerUpdateIS(const LinearController<Scalar, XDimisions, UDimisions, PredictLength> &controller)
{
  std::array<Scalar, PredictLength> biasArraySquaredNorm;

  for (int i = 0; i < PredictLength; ++i)
  {
    biasArraySquaredNorm[i] = controller.deltaBiasArray_[i].squareNorm();
  }
  // integrates using the trapezoidal approximation method
  return trapezoidalIntegration(controller.timeStamp_, biasArraySquaredNorm, 0.0);
}