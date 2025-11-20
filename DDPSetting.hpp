#pragma once
#include "SearchStrategyBase.hpp"
#include "Integration.hpp"
#include "SearchStrategySettings.hpp"

/**
 * This structure contains the settings for the DDP algorithm.
 */
template<typename Scalar>
struct DDPSettings {
  /** Maximum number of iterations of DDP. */
  size_t maxNumIterations_ = 15;
  /** This value determines the termination condition based on the minimum relative changes of the cost. */
  Scalar minRelCost_ = 1e-3;
  /** This value determines the tolerance of constraint's ISE (Integral of Square Error). */
  Scalar constraintTolerance_ = 1e-3;

  /** This value determines the absolute tolerance error for ode solvers. */
  Scalar absTolODE_ = 1e-9;
  /** This value determines the relative tolerance error for ode solvers. */
  Scalar relTolODE_ = 1e-6;
  /** This value determines the maximum number of integration points per a second for ode solvers. */
  int maxNumStepsPerSecond_ = 10000;
  /** The integration time step for Riccati equation which is used for fixed timestep integration scheme. */
  Scalar timeStep_ = 1e-2;
  /** The backward pass integrator type: SLQ uses it for solving Riccati equation and ILQR uses it for discretizing LQ approximation. */
  IntegratorType backwardPassIntegratorType_ = IntegratorType::ODE45;

  /** The initial coefficient of the quadratic penalty function in the merit function. It should be greater than one. */
  Scalar constraintPenaltyInitialValue_ = 2.0;
  /** The rate that the coefficient of the quadratic penalty function in the merit function grows. It should be greater than one. */
  Scalar constraintPenaltyIncreaseRate_ = 2.0;

  /** Use either the optimized control policy (true) or the optimized state-input trajectory (false). */
  bool useFeedbackPolicy_ = false;

  /** The risk sensitivity coefficient for risk aware DDP. */
  Scalar riskSensitiveCoeff_ = 0.0;

  /** Determines the strategy for solving the subproblem. There are two choices line-search strategy and levenberg_marquardt strategy. */
  SearchStrategyType strategy_ = SearchStrategyType::LINE_SEARCH;
  /** The line-search strategy settings. */
  LineSearchSettings<Scalar> lineSearch_;
  // /** The levenberg_marquardt strategy settings. */
  // levenberg_marquardt::Settings levenbergMarquardt_;

};  // end of DDP_Settings
