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

#include "HessianCorrection.hpp"
#include <limits>

/**
 * @brief The DDP strategy enum
 * Enum used in selecting either LINE_SEARCH, LEVENBERG_MARQUARDT, or TRUST_REGION strategies.
 */
enum class SearchStrategyType
{
  LINE_SEARCH,
  LEVENBERG_MARQUARDT
};

/**
 * This structure contains the settings for the Line-Search strategy.
 */
template <typename Scalar>
struct SearchStrategyBaseSettings
{
  /** This value determines to display the log output DDP. */
  bool displayInfo = false;
  /** Printing rollout trajectory for debugging. */
  bool debugPrintRollout = false;
  /** This value determines the termination condition based on the minimum relative changes of the cost. */
  Scalar minRelCost{1e-3};
}; // end of Settings

/**
 * This structure contains the settings for the Line-Search strategy.
 */
template <typename Scalar>
struct LineSearchSettings
{
  /** Minimum step length of line-search strategy. */
  Scalar minStepLength = 0.05;
  /** Maximum step length of line-search strategy. */
  Scalar maxStepLength = 1.0;
  /** Line-search strategy contraction rate. */
  Scalar contractionRate = 0.5;
  /** Armijo coefficient, c defined as f(u + a*p) < f(u) + c*a dfdu.dot(p)  */
  Scalar armijoCoefficient = 1e-4;
  /** The Hessian correction strategy. */
  HessianCorrectionStrategy hessianCorrectionStrategy = HessianCorrectionStrategy::DIAGONAL_SHIFT;
  /** The multiple used for correcting the Hessian for numerical stability of the Riccati backward pass.*/
  Scalar hessianCorrectionMultiple = 1e-6;
}; // end of Settings

/**
 * This structure contains the settings for the Levenberg-Marquardt strategy.
 */
template <typename Scalar>
struct LevenbergMarquardtSettings
{
  /** Minimum pho (the ratio between actual reduction and predicted reduction) to accept the iteration's solution.
   * minAcceptedPho_ should be [0, 0.25);
   * */
  Scalar minAcceptedPho = 0.25;
  /** The default ratio of geometric progression for Riccati multiple. */
  Scalar riccatiMultipleDefaultRatio = 2.0;
  /** The default scalar-factor of geometric progression for Riccati multiple. */
  Scalar riccatiMultipleDefaultFactor = 1e-6;
  /** Maximum number of successive rejections of the iteration's solution. */
  int maxNumSuccessiveRejections = 5;
}; // end of Settings