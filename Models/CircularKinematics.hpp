#pragma once
#include "Cost.hpp"
#include "OptimalControlProblem.hpp"
#include "QuadraticPenalty.hpp"
#include "SystemDynamicsBase.hpp"
/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/**
 * This example defines an optimal control problem where a kinematically modeled
 * particle is supposed to orbit a unite circle (defined as a constraint) with
 * velocity of 1[m/s] (defined as a cost).
 */
namespace circular_model {
static constexpr int STATE_DIM = 2;
static constexpr int INPUT_DIM = 2;
template <typename Scalar>
class CircularKinematicsSystem final
    : public SystemDynamicsBase<Scalar, STATE_DIM, INPUT_DIM> {
 public:
  CircularKinematicsSystem() = default;
  ~CircularKinematicsSystem() override = default;

  Vector<Scalar, STATE_DIM> computeFlowMap(
      Scalar t, const Vector<Scalar, STATE_DIM>& x,
      const Vector<Scalar, INPUT_DIM>& u) const override {
    (void)t;
    (void)x;
    return u;
  }

  VectorFunctionLinearApproximation<Scalar, STATE_DIM, STATE_DIM, INPUT_DIM>
  linearApproximation(Scalar t, const Vector<Scalar, STATE_DIM>& x,
                      const Vector<Scalar, INPUT_DIM>& u) {
    VectorFunctionLinearApproximation<Scalar, STATE_DIM, STATE_DIM, INPUT_DIM>
        dynamics;
    dynamics.f = computeFlowMap(t, x, u);
    dynamics.dfdx.setZero();
    dynamics.dfdu.setIdentity();
    return dynamics;
  }
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/**
 * This example defines an optimal control problem where a kinematically modeled
 * particle is supposed to orbit a unite circle (defined as a constraint) with
 * velocity of 1[m/s] (defined as a cost).
 */
template <typename Scalar, int ArrayLength>
class CircularKinematicsCost
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {
 public:
  CircularKinematicsCost(int cost_number)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number) {
        };

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      const override {
    (void)time;
    (void)timeTrajectory;
    (void)stateTrajectoy;
    (void)inputTrajectory;
    const Scalar angularVelocityError =
        state(0) * input(1) - state(1) * input(0) - Scalar(1.0);
    return Scalar(0.5) * angularVelocityError * angularVelocityError +
           Scalar(0.005) * input.dot(input);
  }
  Scalar getValue(
      int time_index, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      const override {
    (void)time_index;
    (void)timeTrajectory;
    (void)stateTrajectoy;
    (void)inputTrajectory;
    const Scalar angularVelocityError =
        state(0) * input(1) - state(1) * input(0) - Scalar(1.0);
    return Scalar(0.5) * angularVelocityError * angularVelocityError +
           Scalar(0.005) * input.dot(input);
  }
  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      const override {
    (void)time;
    (void)timeTrajectory;
    (void)stateTrajectoy;
    (void)inputTrajectory;
    static_assert(STATE_DIM == 2,
                  "CircularKinematicsCost requires STATE_DIM == 2.");
    static_assert(INPUT_DIM == 2,
                  "CircularKinematicsCost requires INPUT_DIM == 2.");

    ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM> result;

    const Scalar x0 = state(0);
    const Scalar x1 = state(1);
    const Scalar u0 = input(0);
    const Scalar u1 = input(1);

    const Scalar c = x0 * u1 - x1 * u0 - Scalar(1.0);

    result.setZero();
    result.f = Scalar(0.5) * c * c + Scalar(0.005) * input.dot(input);

    result.dfdx(0) = c * u1;
    result.dfdx(1) = -c * u0;

    result.dfdu(0) = -c * x1 + Scalar(0.01) * u0;
    result.dfdu(1) = c * x0 + Scalar(0.01) * u1;

    result.dfdxx(0, 0) = u1 * u1;
    result.dfdxx(0, 1) = -u1 * u0;
    result.dfdxx(1, 0) = -u0 * u1;
    result.dfdxx(1, 1) = u0 * u0;

    result.dfduu(0, 0) = x1 * x1 + Scalar(0.01);
    result.dfduu(0, 1) = -x1 * x0;
    result.dfduu(1, 0) = -x0 * x1;
    result.dfduu(1, 1) = x0 * x0 + Scalar(0.01);

    result.dfdux(0, 0) = -x1 * u1;
    result.dfdux(0, 1) = x1 * u0 - c;
    result.dfdux(1, 0) = x0 * u1 + c;
    result.dfdux(1, 1) = -x0 * u0;
    return result;
  }
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/**
 * This example defines an optimal control problem where a kinematically modeled
 * particle is supposed to orbit a unite circle (defined as a constraint) with
 * velocity of 1[m/s] (defined as a cost).
 */
template <typename Scalar>
class CircularKinematicsConstraints final
    : public StateInputConstraint<Scalar, STATE_DIM, INPUT_DIM> {
 public:
  CircularKinematicsConstraints()
      : StateInputConstraint<Scalar, STATE_DIM, INPUT_DIM>(
            ConstraintOrder::Linear) {}
  ~CircularKinematicsConstraints() override = default;

  Scalar getValue(const Scalar time, const Vector<Scalar, STATE_DIM>& state,
                  const Vector<Scalar, INPUT_DIM>& input) const override {
    (void)time;
    return state.dot(input);
  }

  ScalarFunctionLinearApproximation<Scalar, STATE_DIM, INPUT_DIM>
  getLinearApproximation(
      const Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input) const override {
    ScalarFunctionLinearApproximation<Scalar, STATE_DIM, INPUT_DIM> result;
    result.f = getValue(time, state, input);
    result.dfdx = input;
    result.dfdu = state;
    return result;
  }
};

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar, size_t PredictLength>
inline OptimalControlProblem<Scalar, STATE_DIM, INPUT_DIM, PredictLength, 0, 0,
                             1, 0, 0, 0>&
createCircularKinematicsProblem() {
  using Problem_t = OptimalControlProblem<Scalar, STATE_DIM, INPUT_DIM,
                                          PredictLength, 0, 0, 1, 0, 0, 0>;
  using Cost_t =
      CircularKinematicsCost<Scalar, static_cast<int>(PredictLength + 1)>;
  using Constraint_t = CircularKinematicsConstraints<Scalar>;
  using Penalty_t = QuadraticPenalty<Scalar>;
  using AugmentedLagrangian_t =
      StateInputAugmentedLagrangian<Scalar, STATE_DIM, INPUT_DIM>;

  static CircularKinematicsSystem<Scalar> dynamics;
  static Cost_t cost(0);
  static Constraint_t constraint;
  static Penalty_t penalty(
      typename Penalty_t::Config{Scalar(100.0), Scalar(0.1)});
  static AugmentedLagrangian_t augmentedLagrangian(&constraint, &penalty);
  static Problem_t problem;
  static bool initialized = false;

  if (!initialized) {
    problem.dynamicsPtr = &dynamics;
    problem.cost.add(cost);
    problem.equalityLagrangian.add(&augmentedLagrangian);
    initialized = true;
  }

  return problem;
}

};  // namespace circular_model
