#pragma once

#include "Cost.hpp"
#include "CostCollection.hpp"
#include "Dynamics.hpp"
#include "SystemDynamicsBase.hpp"
#include "Controller.hpp"
#include "Multiplier.hpp"
#include "StateAugmentedLagrangianCollection.hpp"
#include "StateInputAugmentedLagrangianCollection.hpp"

/** Optimal Control Problem definition */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
  int StateEqLagrangianContrainNumbers, int StateInputEqLagrangianContrainNumbers, int StateIneqLagrangianContrainNumbers, int StateInputIneqLagrangianContrainNumbers,
  int FinalStateEqLagrangianContrainNumbers, int FinalStateIneqFinalLagrangianContrainNumbers >
struct OptimalControlProblem
{
  /** Default constructor */
  OptimalControlProblem() = default;

  /** Default destructor */
  ~OptimalControlProblem() = default;

  /** Copy constructor */
  OptimalControlProblem(const OptimalControlProblem& other) = delete;

  /** Copy assignment */
  OptimalControlProblem& operator=(const OptimalControlProblem& rhs) = delete;

  // /** Move constructor */
  // OptimalControlProblem(OptimalControlProblem &&other) noexcept = default;

  // /** Move assignment */
  // OptimalControlProblem &operator=(OptimalControlProblem &&rhs) noexcept = default;

  /* Cost */
  /** Intermediate cost */
  StateInputCostCollection<Scalar, XDimisions, UDimisions, PredictLength + 1> cost;
  /** Intermediate state-only cost */
  StateCostCollection<Scalar, XDimisions, PredictLength + 1> stateCost;

  // /** Final cost */
  StateCostCollection<Scalar, XDimisions, PredictLength + 1> finalCost;

  // target trajectory
  std::array<Scalar, PredictLength + 1> timeTrajectory;
  std::array<Vector<Scalar, XDimisions>, PredictLength + 1> stateTrajectory;
  std::array<Vector<Scalar, UDimisions>, PredictLength + 1> inputTrajectory;

  /* Lagrangians */
  // /** Lagrangian for equality constraints */
  StateInputAugmentedLagrangianCollection<Scalar, XDimisions, UDimisions, StateInputEqLagrangianContrainNumbers> equalityLagrangian;
  // /** Lagrangian for state-only equality constraints */
  StateAugmentedLagrangianCollection<Scalar, XDimisions, StateEqLagrangianContrainNumbers> stateEqualityLagrangian;
  // /** Lagrangian for inequality constraints */
  StateInputAugmentedLagrangianCollection<Scalar, XDimisions, UDimisions, StateInputIneqLagrangianContrainNumbers> inequalityLagrangian;
  // /** Lagrangian for state-only inequality constraints */
  StateAugmentedLagrangianCollection<Scalar, XDimisions, StateIneqLagrangianContrainNumbers> stateInequalityLagrangian;

  /** Lagrangian for final equality constraints */
  StateAugmentedLagrangianCollection<Scalar, XDimisions, FinalStateEqLagrangianContrainNumbers> finalEqualityLagrangian;
  /** Lagrangian for final inequality constraints */
  StateAugmentedLagrangianCollection<Scalar, XDimisions, FinalStateIneqFinalLagrangianContrainNumbers> finalInequalityLagrangian;

  /* Dynamics */
  /** System dynamics pointer */
  SystemDynamicsBase<Scalar, XDimisions, UDimisions>* dynamicsPtr;
};
