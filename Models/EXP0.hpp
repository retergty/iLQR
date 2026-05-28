#pragma once

#include <cstddef>

#include "Cost/QuadraticStateCost.hpp"
#include "Dynamics/LinearSystemDynamics.hpp"
#include "OptimalControl/OptimalControlProblem.hpp"
#include "iLQR/iLQRDescriptor.hpp"

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
            Matrix<Scalar, STATE_DIM, STATE_DIM>{{Scalar(0.6), Scalar(1.2)},
                                                 {Scalar(-0.8), Scalar(3.4)}},
            Matrix<Scalar, STATE_DIM, INPUT_DIM>{Scalar(1.0), Scalar(1.0)}) {}
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
            Matrix<Scalar, STATE_DIM, STATE_DIM>{{Scalar(2.0), Scalar(0.0)},
                                                 {Scalar(0.0), Scalar(2.0)}},
            Matrix<Scalar, INPUT_DIM, INPUT_DIM>::Identity(), 0) {}
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
            Scalar(2) * Matrix<Scalar, STATE_DIM, STATE_DIM>::Identity(), 0) {}
  ~EXP0_FinalCost() override = default;

 private:
  EXP0_FinalCost(const EXP0_FinalCost& other) = default;
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar>
Vector<Scalar, STATE_DIM> getExp0TargetState() {
  return Vector<Scalar, STATE_DIM>{Scalar(4.0), Scalar(2.0)};
}

template <typename Scalar>
Vector<Scalar, INPUT_DIM> getExp0TargetInput() {
  return Vector<Scalar, INPUT_DIM>::Zero();
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar, size_t PredictLength>
inline OptimalControlProblem<
    Scalar,
    TranscriptionConfig<Dimensions<STATE_DIM, INPUT_DIM>,
                        Horizon<PredictLength>>,
    ConstraintConfig<>>&
createExp0Problem() {
  using Transcription_t = TranscriptionConfig<Dimensions<STATE_DIM, INPUT_DIM>,
                                              Horizon<PredictLength>>;
  using Problem_t =
      OptimalControlProblem<Scalar, Transcription_t, ConstraintConfig<>>;
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
