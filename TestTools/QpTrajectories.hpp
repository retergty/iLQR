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

namespace qp_solver {

/** A time, state, input trajectory. The last timepoint has only a state, no
 * input */
template <typename Scalar, int XDim, int UDim, size_t PredictLength>
struct ContinuousTrajectory {
  using TimeTrajectory_t = std::array<Scalar, PredictLength + 1>;
  using StateTrajectory_t = std::array<Vector<Scalar, XDim>, PredictLength + 1>;
  using InputTrajectory_t = std::array<Vector<Scalar, UDim>, PredictLength>;
  /** time trajectory, size N+1 */
  TimeTrajectory_t timeTrajectory;
  /** trajectory of state vectors, size N+1 */
  StateTrajectory_t stateTrajectory;
  /** trajectory of input vectors, size N */
  InputTrajectory_t inputTrajectory;
};

/** Adds state and inputs of two trajectories, time is not added. */
template <typename Scalar, int XDim, int UDim, size_t PredictLength>
ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>
operator+(const ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> &lhs,
          const ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> &rhs) {
  // Copy lhs into sum
  ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> sum(lhs);

  for (size_t k = 0; k < sum.inputTrajectory.size(); ++k) {
    sum.inputTrajectory[k] += rhs.inputTrajectory[k];
  }

  // Sum states
  for (size_t k = 0; k < sum.stateTrajectory.size(); ++k) {
    sum.stateTrajectory[k] += rhs.stateTrajectory[k];
  }
  return sum;
}

/** Reference to a point along a trajectory. Does not own the state-input data.
 */
template <typename Scalar, int XDim, int UDim> struct TrajectoryRef {
  using StateVector_t = Vector<Scalar, XDim>;
  using InputVector_t = Vector<Scalar, UDim>;
  /** time */
  Scalar t;
  /** state */
  const StateVector_t &x;
  /** input */
  const InputVector_t &u;
};

/** Reference to the state at a point along a trajectory. Does not own the state
 * data. */
template <typename Scalar, int XDim> struct StateTrajectoryRef {
  using StateVector_t = Vector<Scalar, XDim>;
  /** time */
  Scalar t;
  /** state */
  const StateVector_t &x;
};

} // namespace qp_solver
