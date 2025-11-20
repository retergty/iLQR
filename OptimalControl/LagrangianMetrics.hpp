#pragma once
#include "Types.hpp"

/** The structure contains a term's constraint vector and its associated penalty */
template<typename Scalar>
struct LagrangianMetrics {
  LagrangianMetrics() : LagrangianMetrics(0, 0) {}
  LagrangianMetrics(Scalar penaltyArg, Scalar constraintArg) : penalty(penaltyArg), constraint(constraintArg) {}

  Scalar penalty;
  Scalar constraint;
};