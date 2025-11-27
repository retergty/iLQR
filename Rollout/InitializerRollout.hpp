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

#include "RolloutBase.hpp"
#include "Initializer.hpp"
/**
 * This class is an interface class for forward rollout of the initializer.
 */
template <typename Scalar, int XDimisions, int UDimisions>
class InitializerRollout : RolloutBase<Scalar, XDimisions, UDimisions>
{
public:
  using Initializer_t = Initializer<Scalar, XDimisions, UDimisions>;
  using RolloutTrajectoryPointer_t = typename RolloutBase<Scalar, XDimisions, UDimisions>::RolloutTrajectoryPointer_t;
  /**
   * Constructor.
   *
   * @param [in] initializer: The initializer for the state and the input.
   * @param [in] rolloutSettings: The rollout settings.
   */
  explicit InitializerRollout(Initializer_t &initializer, const Scalar timeStep) : initializer_(initializer)
  {
    this->rolloutSettings_.timeStep = timeStep;
  }

  ~InitializerRollout() override = default;

  int run(const Scalar initTime, const Vector<Scalar, XDimisions> &initState, const Scalar finalTime, ControllerBase<Scalar, XDimisions, UDimisions> *controller,
          RolloutTrajectoryPointer_t &trajectory) override
  {
    assert(finalTime > initTime);
    (void)controller;

    // Ensure that finalTime is included by adding a fraction of dt such that: N * dt <= finalTime < (N + 1) * dt.
    Scalar finalTimeLocal = finalTime + 0.1 * this->settings().timeStep;
    Scalar t = initTime;
    const Scalar timeStep = this->settings().timeStep;

    Vector<Scalar, XDimisions> state = initState;
    Vector<Scalar, XDimisions> nextState;

    const size_t numSteps = (finalTimeLocal - initTime) / timeStep;

    for (size_t i = 0; i < numSteps + 1; i++)
    {
      assert(i < trajectory.maxLength);
      initializer_.compute(t, state, t + timeStep, trajectory.inputTrajectory[i], nextState);
      trajectory.timeTrajectory[i] = t;
      trajectory.stateTrajectory[i] = state;
      state = nextState;
      t += timeStep;
    } // end of i loop
    return numSteps;
  }

private:
  Initializer_t &initializer_;
};