#pragma once
#include "Cost.hpp"
#include "OptimalControlProblem.hpp"
#include "QuadraticPenalty.hpp"
#include "SmoothAbsolutePenalty.hpp"
#include "StateInputAugmentedLagrangian.hpp"
#include "StateInputConstraint.hpp"
#include "SystemDynamicsBase.hpp"
#include "iLQRDescriptor.hpp"

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/**
 * 该示例定义了一个最优控制问题，其中运动学建模的
 * 粒子需要以 1[m/s] 速度绕单位圆运动（单位圆定义为约束），
 * 速度要求定义为代价。
 */
namespace circular_model {
static constexpr int STATE_DIM = 2;
static constexpr int INPUT_DIM = 2;

enum class AugmentedLagrangianType {
  Quadratic,
  QuadraticStrong,
  SmoothAbsolute
};

template <typename Scalar, size_t PredictLength>
using CircularKinematicsProblem = OptimalControlProblem<
    Scalar,
    TranscriptionConfig<Dimensions<STATE_DIM, INPUT_DIM>,
                        Horizon<PredictLength>>,
    ConstraintConfig<
        StateConstraintConfig<ConstraintLayout<>>,
        StateInputConstraintConfig<ConstraintLayout<
            ConstraintGroupLayout<ConstraintTerm<1>>, ConstraintGroupLayout<>>>,
        FinalStateConstraintConfig<ConstraintLayout<>>>>;

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
 * 该示例定义了一个最优控制问题，其中运动学建模的
 * 粒子需要以 1[m/s] 速度绕单位圆运动（单位圆定义为约束），
 * 速度要求定义为代价。
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
 * 该示例定义了一个最优控制问题，其中运动学建模的
 * 粒子需要以 1[m/s] 速度绕单位圆运动（单位圆定义为约束），
 * 速度要求定义为代价。
 */
template <typename Scalar>
class CircularKinematicsConstraints final
    : public StateInputConstraint<Scalar, STATE_DIM, INPUT_DIM, 1> {
 public:
  CircularKinematicsConstraints()
      : StateInputConstraint<Scalar, STATE_DIM, INPUT_DIM, 1>(
            ConstraintOrder::Linear) {}
  ~CircularKinematicsConstraints() override = default;

  Vector<Scalar, 1> getValue(
      const Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input) const override {
    (void)time;
    Vector<Scalar, 1> value;
    value(0) = state.dot(input);
    return value;
  }

  VectorFunctionLinearApproximation<Scalar, 1, STATE_DIM, INPUT_DIM>
  getLinearApproximation(
      const Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input) const override {
    VectorFunctionLinearApproximation<Scalar, 1, STATE_DIM, INPUT_DIM> result;
    result.f = getValue(time, state, input);
    result.dfdx = input.transpose();
    result.dfdu = state.transpose();
    return result;
  }
};

template <typename Scalar, size_t PredictLength, AugmentedLagrangianType Type,
          typename Penalty_t>
inline CircularKinematicsProblem<Scalar, PredictLength>&
createCircularKinematicsProblemWithPenalty(Penalty_t& penalty) {
  using Cost_t =
      CircularKinematicsCost<Scalar, static_cast<int>(PredictLength + 1)>;
  using Constraint_t = CircularKinematicsConstraints<Scalar>;
  using AugmentedLagrangian_t =
      StateInputAugmentedLagrangian<Scalar, STATE_DIM, INPUT_DIM, 1>;
  using Problem_t = CircularKinematicsProblem<Scalar, PredictLength>;

  static CircularKinematicsSystem<Scalar> dynamics;
  static Cost_t cost(0);
  static Constraint_t constraint;
  static AugmentedLagrangian_t augmentedLagrangian(&constraint, &penalty);
  static Problem_t problem;
  static bool initialized = false;

  if (!initialized) {
    problem.dynamicsPtr = &dynamics;
    problem.cost.add(cost);
    problem.equalityLagrangian.template set<0>(&augmentedLagrangian);
    initialized = true;
  }

  return problem;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
template <typename Scalar, size_t PredictLength>
inline CircularKinematicsProblem<Scalar, PredictLength>&
createCircularKinematicsProblemImpl(AugmentedLagrangianType type) {
  using QuadraticPenalty_t = QuadraticPenalty<Scalar>;
  using SmoothAbsolutePenalty_t = SmoothAbsolutePenalty<Scalar>;

  switch (type) {
    case AugmentedLagrangianType::Quadratic: {
      static QuadraticPenalty_t penalty(
          typename QuadraticPenalty_t::Config{Scalar(100.0), Scalar(0.1)});
      return createCircularKinematicsProblemWithPenalty<
          Scalar, PredictLength, AugmentedLagrangianType::Quadratic>(penalty);
    }
    case AugmentedLagrangianType::QuadraticStrong: {
      static QuadraticPenalty_t penalty(
          typename QuadraticPenalty_t::Config{Scalar(150.0), Scalar(0.1)});
      return createCircularKinematicsProblemWithPenalty<
          Scalar, PredictLength, AugmentedLagrangianType::QuadraticStrong>(
          penalty);
    }
    case AugmentedLagrangianType::SmoothAbsolute: {
      static SmoothAbsolutePenalty_t penalty(
          typename SmoothAbsolutePenalty_t::Config{Scalar(0.1), Scalar(0.1),
                                                   Scalar(0.1)});
      return createCircularKinematicsProblemWithPenalty<
          Scalar, PredictLength, AugmentedLagrangianType::SmoothAbsolute>(
          penalty);
    }
  }

  static QuadraticPenalty_t fallbackPenalty(
      typename QuadraticPenalty_t::Config{Scalar(100.0), Scalar(0.1)});
  return createCircularKinematicsProblemWithPenalty<
      Scalar, PredictLength, AugmentedLagrangianType::Quadratic>(
      fallbackPenalty);
}

template <typename Scalar, size_t PredictLength>
inline CircularKinematicsProblem<Scalar, PredictLength>&
createCircularKinematicsProblem(
    AugmentedLagrangianType type = AugmentedLagrangianType::Quadratic) {
  return createCircularKinematicsProblemImpl<Scalar, PredictLength>(type);
}

};  // namespace circular_model
