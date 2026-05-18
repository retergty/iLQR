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

#include <cstddef>

#include "LinearSystemDynamics.hpp"
#include "OptimalControlProblem.hpp"
#include "QuadraticStateCost.hpp"

namespace exp0 {

static constexpr int STATE_DIM = 2;
static constexpr int INPUT_DIM = 1;

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar>
class EXP0_Sys1 final
    : public LinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM> {
 public:
  EXP0_Sys1()
      : LinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM>(
            (Matrix<Scalar, STATE_DIM, STATE_DIM>() << Scalar(0.6), Scalar(1.2),
             Scalar(-0.8), Scalar(3.4))
                .finished(),
            (Matrix<Scalar, STATE_DIM, INPUT_DIM>() << Scalar(1.0), Scalar(1.0))
                .finished()) {}
  ~EXP0_Sys1() override = default;

 private:
  EXP0_Sys1(const EXP0_Sys1& other) = default;
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar, int ArrayLength>
class EXP0_Cost final : public QuadraticStateInputCost<Scalar, STATE_DIM,
                                                       INPUT_DIM, ArrayLength> {
 public:
  EXP0_Cost()
      : QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(
            (Matrix<Scalar, STATE_DIM, STATE_DIM>() << Scalar(0.0), Scalar(0.0),
             Scalar(0.0), Scalar(1.0))
                .finished(),
            Matrix<Scalar, INPUT_DIM, INPUT_DIM>::Identity()) {}
  ~EXP0_Cost() override = default;

 private:
  EXP0_Cost(const EXP0_Cost& other) = default;
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar, int ArrayLength>
class EXP0_FinalCost final
    : public QuadraticStateCost<Scalar, STATE_DIM, ArrayLength> {
 public:
  EXP0_FinalCost()
      : QuadraticStateCost<Scalar, STATE_DIM, ArrayLength>(
            Matrix<Scalar, STATE_DIM, STATE_DIM>::Identity()) {}
  ~EXP0_FinalCost() override = default;

 private:
  EXP0_FinalCost(const EXP0_FinalCost& other) = default;
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar>
Vector<Scalar, STATE_DIM> getExp0TargetState() {
  return (Vector<Scalar, STATE_DIM>() << Scalar(4.0), Scalar(2.0)).finished();
}

template <typename Scalar>
Vector<Scalar, INPUT_DIM> getExp0TargetInput() {
  return Vector<Scalar, INPUT_DIM>::Zero();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar, size_t PredictLength>
inline OptimalControlProblem<Scalar, STATE_DIM, INPUT_DIM, PredictLength, 0, 0,
                             0, 0, 0, 0>&
createExp0Problem() {
  using Problem_t = OptimalControlProblem<Scalar, STATE_DIM, INPUT_DIM,
                                          PredictLength, 0, 0, 0, 0, 0, 0>;
  using Cost_t = EXP0_Cost<Scalar, static_cast<int>(PredictLength + 1)>;
  using FinalCost_t =
      EXP0_FinalCost<Scalar, static_cast<int>(PredictLength + 1)>;

  static EXP0_Sys1<Scalar> dynamics;
  static Cost_t cost;
  static FinalCost_t finalCost;
  static Problem_t problem;
  static bool initialized = false;

  if (!initialized) {
    problem.dynamicsPtr = &dynamics;
    problem.cost.add(cost);
    problem.finalCost.add(finalCost);
    initialized = true;
  }
  return problem;
}

}  // namespace exp0
