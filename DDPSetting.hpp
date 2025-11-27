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
  size_t maxNumIterations_ = 10;
  /** This value determines the termination condition based on the minimum relative changes of the cost. */
  Scalar minRelCost_ = 1e-3;
  /** This value determines the tolerance of constraint's ISE (Integral of Square Error). */
  Scalar constraintTolerance_ = 1e-3;

  /** The integration time step for Riccati equation which is used for fixed timestep integration scheme. */
  Scalar timeStep_ = 1e-2;

  /** Use either the optimized control policy (true) or the optimized state-input trajectory (false). */
  bool useFeedbackPolicy_ = false;

  /** The risk sensitivity coefficient for risk aware DDP. */
  Scalar riskSensitiveCoeff_ = 0.0;

  /** Determines the strategy for solving the subproblem. There are two choices line-search strategy and levenberg_marquardt strategy. */
  SearchStrategyType strategy_ = SearchStrategyType::LINE_SEARCH;
  /** The line-search strategy settings. */
  LineSearchSettings<Scalar> lineSearch_{};

};  // end of DDP_Settings
