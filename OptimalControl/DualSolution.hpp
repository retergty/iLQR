#pragma once
#include <array>
#include "Multiplier.hpp"
#include "LinearInterpolation.hpp"

template <typename Scalar, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains, int FinalStateEqConstrains, int FinalStateIneqConstrains, size_t PredictLength>
struct DualSolution
{
  using IntermediateMultiplierCollection_t = MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>;
  using FinalMultiplierCollection_t = MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0>;

  std::array<Scalar, PredictLength> timeTrajectory;

  FinalMultiplierCollection_t final;
  std::array<IntermediateMultiplierCollection_t, PredictLength> intermediates;

  /** Exchanges the content of DualSolution */
  void swap(DualSolution &other)
  {
    timeTrajectory.swap(other.timeTrajectory);
    intermediates.swap(other.intermediates);
    final.swap(other.final);
  }
};

/**
 * Note that this view of DualSolution does not have access to its timestamp, signaling that functions
 * operating on it cannot modify its timestamp.
 */
template <typename Scalar, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains, int FinalStateEqConstrains, int FinalStateIneqConstrains, size_t PredictLength>
struct DualSolutionRef
{
  using DualSolution_t = DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength>;
  using IntermediateMultiplierCollection_t = MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>;
  using FinalMultiplierCollection_t = MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0>;

  DualSolutionRef(DualSolution_t &dualSolution) : DualSolutionRef(dualSolution.final, dualSolution.intermediates) {}

  DualSolutionRef(FinalMultiplierCollection_t &finalRef, std::array<IntermediateMultiplierCollection_t, PredictLength> &intermediatesRef)
      : final(finalRef), intermediates(intermediatesRef) {}

  FinalMultiplierCollection_t &final;
  std::array<IntermediateMultiplierCollection_t, PredictLength> &intermediates;
};

/**
 * Calculates the intermediate dual solution at the given time.
 *
 * @param [in] dualSolution: The dual solution
 * @param [in] time: The inquiry time
 * @return The collection of multipliers associated to state/state-input, equality/inequality Lagrangian terms.
 */
template <typename Scalar, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains, int FinalStateEqConstrains, int FinalStateIneqConstrains, size_t PredictLength>
inline MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>
getIntermediateDualSolutionAtTime(const DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength> &dualSolution,
                                  Scalar time)
{
  const std::pair<int, Scalar> indexAlpha = LinearInterpolation::timeSegment(time, dualSolution.timeTrajectory);
  return LinearInterpolation::interpolate(indexAlpha, dualSolution.intermediates);
}