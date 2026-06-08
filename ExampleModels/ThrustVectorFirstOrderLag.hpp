#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <limits>

#include "Cost/Cost.hpp"
#include "Dynamics/DiscreteSystemDynamicsBase.hpp"
#include "Initialization/Initializer.hpp"
#include "Misc/LinearInterpolation.hpp"
#include "OptimalControl/OptimalControlProblem.hpp"
#include "iLQR/DDPSetting.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"
#include "iLQR/iLQR.hpp"
#include "iLQR/iLQRDescriptor.hpp"

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/**
 * 一阶低通执行器版本的矢量推力无人机离散动力学。
 *
 * 状态为 9 维：
 *   [NED 速度, NED 实际生效总加速度, 上一拍 NED 命令总加速度]
 * 输入为 3 维：
 *   当前 NED 命令总加速度增量
 */

namespace thrust_vector_lag {
static constexpr int STATE_DIM = 9;
static constexpr int INPUT_DIM = 3;
template <typename Scalar, size_t PredictLength>
using ThrustVectorProblem = OptimalControlProblem<
    Scalar,
    TranscriptionConfig<Dimensions<STATE_DIM, INPUT_DIM>,
                        Horizon<PredictLength>, DiscreteDynamics>,
    ConstraintConfig<>>;

template <typename Scalar>
class ThrustVectorDynamicSystem final
    : public DiscreteSystemDynamicsBase<Scalar, STATE_DIM, INPUT_DIM> {
 public:
  using LinearApproximation_t =
      typename DiscreteSystemDynamicsBase<Scalar, STATE_DIM,
                                          INPUT_DIM>::LinearApproximation_t;

  struct PreCompCache {
    PreCompCache() {
      approximation.setZero();
      dt = 0;
      alpha = 1;
      oneMinusAlpha = 0;
      velocityEffectiveCoefficient = 0;
      velocityCommandCoefficient = 0;
      dirty = true;
    }

    LinearApproximation_t approximation;
    Scalar dt;
    Scalar alpha;
    Scalar oneMinusAlpha;
    Scalar velocityEffectiveCoefficient;
    Scalar velocityCommandCoefficient;
    bool dirty;
  };

  explicit ThrustVectorDynamicSystem(const Scalar alpha = Scalar(1))
      : alpha_(alpha) {
    assert(alpha_ >= Scalar(0));
    assert(alpha_ <= Scalar(1));
  }

  ~ThrustVectorDynamicSystem() override = default;

  void updateCache(const Scalar dt) {
    if (cache_.dirty ||
        std::abs(cache_.dt - dt) > std::numeric_limits<Scalar>::epsilon() ||
        std::abs(cache_.alpha - alpha_) >
            std::numeric_limits<Scalar>::epsilon()) {
      const Scalar alpha = alpha_;
      const Scalar oneMinusAlpha = Scalar(1) - alpha;
      const Scalar velocityEffectiveCoefficient = dt * oneMinusAlpha;
      const Scalar velocityCommandCoefficient = dt * alpha;

      cache_.approximation.setZero();

      for (int i = 0; i < 3; ++i) {
        cache_.approximation.dfdx(i, i) = Scalar(1);
        cache_.approximation.dfdx(i, i + 3) = velocityEffectiveCoefficient;
        cache_.approximation.dfdx(i, i + 6) = velocityCommandCoefficient;

        cache_.approximation.dfdx(i + 3, i + 3) = oneMinusAlpha;
        cache_.approximation.dfdx(i + 3, i + 6) = alpha;

        cache_.approximation.dfdx(i + 6, i + 6) = Scalar(1);

        cache_.approximation.dfdu(i, i) = velocityCommandCoefficient;
        cache_.approximation.dfdu(i + 3, i) = alpha;
        cache_.approximation.dfdu(i + 6, i) = Scalar(1);
      }

      cache_.dt = dt;
      cache_.alpha = alpha;
      cache_.oneMinusAlpha = oneMinusAlpha;
      cache_.velocityEffectiveCoefficient = velocityEffectiveCoefficient;
      cache_.velocityCommandCoefficient = velocityCommandCoefficient;
      cache_.dirty = false;
    }
  }

  LinearApproximation_t linearApproximation(Scalar t,
                                            const Vector<Scalar, STATE_DIM>& x,
                                            const Vector<Scalar, INPUT_DIM>& u,
                                            Scalar dt) override {
    (void)t;

    updateCache(dt);
    LinearApproximation_t approximation = cache_.approximation;
    approximation.f = computeMapCached(x, u);
    return approximation;
  }

  LinearApproximation_t deviationLinearApproximation(
      Scalar t, const Vector<Scalar, STATE_DIM>& x,
      const Vector<Scalar, INPUT_DIM>& u, Scalar dt) override {
    (void)t;
    (void)x;
    (void)u;

    updateCache(dt);
    return cache_.approximation;
  }

  Vector<Scalar, STATE_DIM> computeMap(Scalar t,
                                       const Vector<Scalar, STATE_DIM>& x,
                                       const Vector<Scalar, INPUT_DIM>& u,
                                       Scalar dt) override {
    (void)t;

    updateCache(dt);
    return computeMapCached(x, u);
  }

 private:
  Vector<Scalar, STATE_DIM> computeMapCached(
      const Vector<Scalar, STATE_DIM>& x,
      const Vector<Scalar, INPUT_DIM>& u) const {
    Vector<Scalar, STATE_DIM> nextState;
    for (int i = 0; i < 3; ++i) {
      const Scalar commandAcceleration = x(i + 6) + u(i);
      nextState(i) = x(i) + cache_.velocityEffectiveCoefficient * x(i + 3) +
                     cache_.velocityCommandCoefficient * commandAcceleration;
      nextState(i + 3) =
          cache_.oneMinusAlpha * x(i + 3) + cache_.alpha * commandAcceleration;
      nextState(i + 6) = commandAcceleration;
    }
    return nextState;
  }

  Scalar alpha_;
  PreCompCache cache_;
};

template <typename Scalar, int ArrayLength>
class ThrustVectorTrackCost final
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {
 public:
  ThrustVectorTrackCost(const Matrix<Scalar, STATE_DIM, STATE_DIM>& Q,
                        const Matrix<Scalar, INPUT_DIM, INPUT_DIM>& R,
                        int cost_number)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number),
        Qv_(Q.template topLeftCorner<3, 3>()),
        R_(R) {
    approximation_.setZero();
    approximation_.dfdxx.template topLeftCorner<3, 3>() = Qv_;
    approximation_.dfduu = R_;
  }

  ~ThrustVectorTrackCost() override = default;

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    const Vector<Scalar, 3> velocityDeviation =
        state.template head<3>() -
        interpolateVelocityReference(indexAlpha, stateTrajectory);
    const Vector<Scalar, INPUT_DIM> inputDeviation =
        input - interpolateInputReference(indexAlpha, inputTrajectory);
    const Vector<Scalar, 3> weightedVelocityDeviation = Qv_ * velocityDeviation;
    const Vector<Scalar, INPUT_DIM> weightedInputDeviation =
        R_ * inputDeviation;

    return Scalar(0.5) * velocityDeviation.dot(weightedVelocityDeviation) +
           Scalar(0.5) * inputDeviation.dot(weightedInputDeviation);
  }

  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    const Vector<Scalar, 3> velocityDeviation =
        state.template head<3>() -
        interpolateVelocityReference(indexAlpha, stateTrajectory);
    const Vector<Scalar, INPUT_DIM> inputDeviation =
        input - interpolateInputReference(indexAlpha, inputTrajectory);
    const Vector<Scalar, 3> weightedVelocityDeviation = Qv_ * velocityDeviation;
    const Vector<Scalar, INPUT_DIM> weightedInputDeviation =
        R_ * inputDeviation;

    approximation_.f =
        Scalar(0.5) * velocityDeviation.dot(weightedVelocityDeviation) +
        Scalar(0.5) * inputDeviation.dot(weightedInputDeviation);
    approximation_.dfdx.template head<3>() = weightedVelocityDeviation;
    approximation_.dfdu = weightedInputDeviation;
    return approximation_;
  }

 private:
  ThrustVectorTrackCost(const ThrustVectorTrackCost& other) = default;

  Vector<Scalar, 3> interpolateVelocityReference(
      const std::pair<int, Scalar>& indexAlpha,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory)
      const {
    return LinearInterpolation::interpolate(
        indexAlpha, stateTrajectory,
        [](const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& trajectory,
           size_t index) -> Vector<Scalar, 3> {
          return trajectory[index].template head<3>();
        });
  }

  Vector<Scalar, INPUT_DIM> interpolateInputReference(
      const std::pair<int, Scalar>& indexAlpha,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      const {
    return LinearInterpolation::interpolate(indexAlpha, inputTrajectory);
  }

  Matrix<Scalar, 3, 3> Qv_;
  Matrix<Scalar, INPUT_DIM, INPUT_DIM> R_;
  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
      approximation_;
};

template <typename Scalar, int ArrayLength>
class ThrustVectorDiagonalTrackCost final
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {
 public:
  ThrustVectorDiagonalTrackCost(const Matrix<Scalar, STATE_DIM, STATE_DIM>& Q,
                                const Matrix<Scalar, INPUT_DIM, INPUT_DIM>& R,
                                int cost_number)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number) {
    approximation_.setZero();
    for (int i = 0; i < 3; ++i) {
      QvDiagonal_(i) = Q(i, i);
      RDiagonal_(i) = R(i, i);
      approximation_.dfdxx(i, i) = QvDiagonal_(i);
      approximation_.dfduu(i, i) = RDiagonal_(i);
    }
  }

  ~ThrustVectorDiagonalTrackCost() override = default;

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    return computeDiagonalCost<false>(indexAlpha, state, input, stateTrajectory,
                                      inputTrajectory);
  }

  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    approximation_.f = computeDiagonalCost<true>(
        indexAlpha, state, input, stateTrajectory, inputTrajectory);
    return approximation_;
  }

 private:
  ThrustVectorDiagonalTrackCost(const ThrustVectorDiagonalTrackCost& other) =
      default;

  Scalar interpolateStateReferenceComponent(
      const std::pair<int, Scalar>& indexAlpha,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      int dim) const {
    if constexpr (ArrayLength > 1) {
      const int index = indexAlpha.first;
      const Scalar alpha = indexAlpha.second;
      return alpha * stateTrajectory[index](dim) +
             (Scalar(1) - alpha) * stateTrajectory[index + 1](dim);
    } else {
      return stateTrajectory[0](dim);
    }
  }

  Scalar interpolateInputReferenceComponent(
      const std::pair<int, Scalar>& indexAlpha,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory,
      int dim) const {
    if constexpr (ArrayLength > 1) {
      const int index = indexAlpha.first;
      const Scalar alpha = indexAlpha.second;
      return alpha * inputTrajectory[index](dim) +
             (Scalar(1) - alpha) * inputTrajectory[index + 1](dim);
    } else {
      return inputTrajectory[0](dim);
    }
  }

  template <bool UpdateApproximation>
  Scalar computeDiagonalCost(
      const std::pair<int, Scalar>& indexAlpha,
      const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>&
          inputTrajectory) {
    Scalar value = 0;
    for (int i = 0; i < 3; ++i) {
      const Scalar velocityDeviation =
          state(i) -
          interpolateStateReferenceComponent(indexAlpha, stateTrajectory, i);
      const Scalar inputDeviation =
          input(i) -
          interpolateInputReferenceComponent(indexAlpha, inputTrajectory, i);

      value += QvDiagonal_(i) * velocityDeviation * velocityDeviation +
               RDiagonal_(i) * inputDeviation * inputDeviation;

      if constexpr (UpdateApproximation) {
        approximation_.dfdx(i) = QvDiagonal_(i) * velocityDeviation;
        approximation_.dfdu(i) = RDiagonal_(i) * inputDeviation;
      }
    }
    return Scalar(0.5) * value;
  }

  Vector<Scalar, 3> QvDiagonal_;
  Vector<Scalar, INPUT_DIM> RDiagonal_;
  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
      approximation_;
};

template <typename Scalar, int ArrayLength>
class ThrustVectorTrackFinalCost final
    : public StateCost<Scalar, STATE_DIM, ArrayLength> {
 public:
  ThrustVectorTrackFinalCost(const Matrix<Scalar, STATE_DIM, STATE_DIM>& Q,
                             int cost_number)
      : StateCost<Scalar, STATE_DIM, ArrayLength>(cost_number),
        Qv_(Q.template topLeftCorner<3, 3>()) {
    approximation_.setZero();
    approximation_.dfdxx.template topLeftCorner<3, 3>() = Qv_;
  }

  ~ThrustVectorTrackFinalCost() override = default;

  Scalar getValue(Scalar time, const Vector<Scalar, STATE_DIM>& state,
                  const std::array<Scalar, ArrayLength>& timeTrajectory,
                  const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>&
                      stateTrajectory) override {
    const Vector<Scalar, 3> velocityDeviation =
        state.template head<3>() -
        interpolateVelocityReference(time, timeTrajectory, stateTrajectory);
    return Scalar(0.5) * velocityDeviation.dot(Qv_ * velocityDeviation);
  }

  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, 0>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory)
      override {
    const Vector<Scalar, 3> velocityDeviation =
        state.template head<3>() -
        interpolateVelocityReference(time, timeTrajectory, stateTrajectory);

    approximation_.f =
        Scalar(0.5) * velocityDeviation.dot(Qv_ * velocityDeviation);
    approximation_.dfdx.template head<3>() = Qv_ * velocityDeviation;
    return approximation_;
  }

 private:
  ThrustVectorTrackFinalCost(const ThrustVectorTrackFinalCost& other) = default;

  Vector<Scalar, 3> interpolateVelocityReference(
      Scalar time, const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory)
      const {
    return LinearInterpolation::interpolate(
        time, timeTrajectory, stateTrajectory,
        [](const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& trajectory,
           size_t index) -> Vector<Scalar, 3> {
          return trajectory[index].template head<3>();
        });
  }

  Matrix<Scalar, 3, 3> Qv_;
  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, 0> approximation_;
};

template <typename Scalar, int ArrayLength>
class ThrustVectorDiagonalTrackFinalCost final
    : public StateCost<Scalar, STATE_DIM, ArrayLength> {
 public:
  ThrustVectorDiagonalTrackFinalCost(
      const Matrix<Scalar, STATE_DIM, STATE_DIM>& Q, int cost_number)
      : StateCost<Scalar, STATE_DIM, ArrayLength>(cost_number) {
    approximation_.setZero();
    for (int i = 0; i < 3; ++i) {
      QvDiagonal_(i) = Q(i, i);
      approximation_.dfdxx(i, i) = QvDiagonal_(i);
    }
  }

  ~ThrustVectorDiagonalTrackFinalCost() override = default;

  Scalar getValue(Scalar time, const Vector<Scalar, STATE_DIM>& state,
                  const std::array<Scalar, ArrayLength>& timeTrajectory,
                  const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>&
                      stateTrajectory) override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    return computeDiagonalCost<false>(indexAlpha, state, stateTrajectory);
  }

  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, 0>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    approximation_.f =
        computeDiagonalCost<true>(indexAlpha, state, stateTrajectory);
    return approximation_;
  }

 private:
  ThrustVectorDiagonalTrackFinalCost(
      const ThrustVectorDiagonalTrackFinalCost& other) = default;

  Scalar interpolateStateReferenceComponent(
      const std::pair<int, Scalar>& indexAlpha,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      int dim) const {
    if constexpr (ArrayLength > 1) {
      const int index = indexAlpha.first;
      const Scalar alpha = indexAlpha.second;
      return alpha * stateTrajectory[index](dim) +
             (Scalar(1) - alpha) * stateTrajectory[index + 1](dim);
    } else {
      return stateTrajectory[0](dim);
    }
  }

  template <bool UpdateApproximation>
  Scalar computeDiagonalCost(const std::pair<int, Scalar>& indexAlpha,
                             const Vector<Scalar, STATE_DIM>& state,
                             const std::array<Vector<Scalar, STATE_DIM>,
                                              ArrayLength>& stateTrajectory) {
    Scalar value = 0;
    for (int i = 0; i < 3; ++i) {
      const Scalar velocityDeviation =
          state(i) -
          interpolateStateReferenceComponent(indexAlpha, stateTrajectory, i);
      value += QvDiagonal_(i) * velocityDeviation * velocityDeviation;

      if constexpr (UpdateApproximation) {
        approximation_.dfdx(i) = QvDiagonal_(i) * velocityDeviation;
      }
    }
    return Scalar(0.5) * value;
  }

  Vector<Scalar, 3> QvDiagonal_;
  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, 0> approximation_;
};

template <typename Scalar, int ArrayLength>
class ThrustDirectionChangeCost final
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {
 public:
  static constexpr Scalar epsilon = Scalar(1e-4);
  static constexpr Scalar Weight = 1;
  static constexpr Scalar MinThrustAccelerationForDirection = Scalar(1e-2);
  static constexpr Scalar GravityAcceleration = Scalar(9.8);

  explicit ThrustDirectionChangeCost(int cost_number = 0)
      : ThrustDirectionChangeCost(Scalar(1), Weight, cost_number) {}

  ThrustDirectionChangeCost(const Scalar alpha, int cost_number)
      : ThrustDirectionChangeCost(alpha, Weight, cost_number) {}

  ThrustDirectionChangeCost(const Scalar alpha, const Scalar weight,
                            int cost_number)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number),
        alpha_(alpha),
        oneMinusAlpha_(Scalar(1) - alpha),
        weight_(weight),
        gravityVector_{Scalar(0), Scalar(0), GravityAcceleration} {
    assert(alpha_ >= Scalar(0));
    assert(alpha_ <= Scalar(1));
    approximation_.setZero();
  }

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    (void)time;
    (void)timeTrajectory;
    (void)inputTrajectory;
    (void)stateTrajectory;

    const Vector<Scalar, 3> lastThrustAcceleration =
        effectiveAcceleration(state) - gravityVector_;
    const Vector<Scalar, 3> thrustAcceleration =
        nextEffectiveAcceleration(state, input) - gravityVector_;
    const Vector<Scalar, 3> lastThrustDirection =
        lastThrustAcceleration /
        std::sqrt(lastThrustAcceleration.dot(lastThrustAcceleration) + epsilon);
    const Vector<Scalar, 3> thrustDirection =
        thrustAcceleration /
        std::sqrt(thrustAcceleration.dot(thrustAcceleration) + epsilon);
    const Vector<Scalar, 3> rk = thrustDirection - lastThrustDirection;
    const Scalar gate =
        lowAccelerationGate(lastThrustAcceleration, thrustAcceleration);
    return Scalar(0.5) * gate * weight_ * rk.dot(rk);
  }

  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    (void)time;
    (void)timeTrajectory;
    (void)inputTrajectory;
    (void)stateTrajectory;

    const Vector<Scalar, 3> lastThrustAcceleration =
        effectiveAcceleration(state) - gravityVector_;
    const Vector<Scalar, 3> thrustAcceleration =
        nextEffectiveAcceleration(state, input) - gravityVector_;
    const Scalar lastThrustAccelerationNorm =
        std::sqrt(lastThrustAcceleration.dot(lastThrustAcceleration) + epsilon);
    const Scalar thrustAccelerationNorm =
        std::sqrt(thrustAcceleration.dot(thrustAcceleration) + epsilon);
    const Vector<Scalar, 3> lastThrustDirection =
        lastThrustAcceleration / lastThrustAccelerationNorm;
    const Vector<Scalar, 3> thrustDirection =
        thrustAcceleration / thrustAccelerationNorm;
    const Vector<Scalar, 3> rk = thrustDirection - lastThrustDirection;

    const Matrix<Scalar, 3, 3> lastThrustAccelerationJacobian =
        computeDirectionJacobian(lastThrustAcceleration,
                                 lastThrustAccelerationNorm);
    const Matrix<Scalar, 3, 3> thrustAccelerationJacobian =
        computeDirectionJacobian(thrustAcceleration, thrustAccelerationNorm);
    const Matrix<Scalar, 3, 3> effectiveAccelerationJacobian =
        oneMinusAlpha_ * thrustAccelerationJacobian -
        lastThrustAccelerationJacobian;
    const Matrix<Scalar, 3, 3> commandAccelerationJacobian =
        alpha_ * thrustAccelerationJacobian;
    const Matrix<Scalar, 3, 3> inputAccelerationJacobian =
        commandAccelerationJacobian;

    const Scalar gate =
        lowAccelerationGate(lastThrustAcceleration, thrustAcceleration);
    const Scalar effectiveWeight = gate * weight_;

    approximation_.f = Scalar(0.5) * effectiveWeight * rk.dot(rk);
    approximation_.dfdx.template segment<3>(3) =
        effectiveWeight * effectiveAccelerationJacobian.transpose() * rk;
    approximation_.dfdx.template segment<3>(6) =
        effectiveWeight * commandAccelerationJacobian.transpose() * rk;
    approximation_.dfdu =
        effectiveWeight * inputAccelerationJacobian.transpose() * rk;

    approximation_.dfdxx.template slice<3, 3>(3, 3) =
        effectiveWeight * effectiveAccelerationJacobian.transpose() *
        effectiveAccelerationJacobian;
    approximation_.dfdxx.template slice<3, 3>(3, 6) =
        effectiveWeight * effectiveAccelerationJacobian.transpose() *
        commandAccelerationJacobian;
    approximation_.dfdxx.template slice<3, 3>(6, 3) =
        effectiveWeight * commandAccelerationJacobian.transpose() *
        effectiveAccelerationJacobian;
    approximation_.dfdxx.template slice<3, 3>(6, 6) =
        effectiveWeight * commandAccelerationJacobian.transpose() *
        commandAccelerationJacobian;
    approximation_.dfduu = effectiveWeight *
                           inputAccelerationJacobian.transpose() *
                           inputAccelerationJacobian;
    approximation_.dfdux.template slice<3, 3>(0, 3) =
        effectiveWeight * inputAccelerationJacobian.transpose() *
        effectiveAccelerationJacobian;
    approximation_.dfdux.template slice<3, 3>(0, 6) =
        effectiveWeight * inputAccelerationJacobian.transpose() *
        commandAccelerationJacobian;

    return approximation_;
  }

 private:
  static Matrix<Scalar, 3, 3> computeDirectionJacobian(
      const Vector<Scalar, 3>& acceleration, const Scalar norm) {
    Matrix<Scalar, 3, 3> jacobian;
    const Scalar invNorm = Scalar(1) / norm;
    const Scalar invNormCubed = invNorm * invNorm * invNorm;

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        jacobian(i, j) = (i == j ? invNorm : Scalar(0)) -
                         acceleration(i) * acceleration(j) * invNormCubed;
      }
    }
    return jacobian;
  }

  static Scalar accelerationGateSigma(const Vector<Scalar, 3>& acceleration) {
    const Scalar accelerationNormSquared = acceleration.squaredNorm();
    const Scalar minAccelerationSquared =
        MinThrustAccelerationForDirection * MinThrustAccelerationForDirection;
    return accelerationNormSquared /
           (accelerationNormSquared + minAccelerationSquared);
  }

  static Scalar lowAccelerationGate(
      const Vector<Scalar, 3>& lastAcceleration,
      const Vector<Scalar, 3>& currentAcceleration) {
    return accelerationGateSigma(lastAcceleration) *
           accelerationGateSigma(currentAcceleration);
  }

  static Vector<Scalar, 3> effectiveAcceleration(
      const Vector<Scalar, STATE_DIM>& state) {
    Vector<Scalar, 3> acceleration;
    for (int i = 0; i < 3; ++i) {
      acceleration(i) = state(i + 3);
    }
    return acceleration;
  }

  Vector<Scalar, 3> nextEffectiveAcceleration(
      const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input) const {
    Vector<Scalar, 3> acceleration;
    for (int i = 0; i < 3; ++i) {
      const Scalar commandAcceleration = state(i + 6) + input(i);
      acceleration(i) =
          oneMinusAlpha_ * state(i + 3) + alpha_ * commandAcceleration;
    }
    return acceleration;
  }

  Scalar alpha_;
  Scalar oneMinusAlpha_;
  const Scalar weight_;
  const Vector<Scalar, 3> gravityVector_;
  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
      approximation_;
};

template <typename Scalar>
struct ThrustVectorOCPSettings {
  // 权重
  Matrix<Scalar, STATE_DIM, STATE_DIM> Q;
  Matrix<Scalar, INPUT_DIM, INPUT_DIM> R;
  Matrix<Scalar, STATE_DIM, STATE_DIM> Qf;
  Scalar weight;

  // 离散一阶低通系数
  Scalar alpha;
};

template <typename Scalar>
struct ThrustVectorILQRSettings {
  DDPSettings<Scalar> ddpSettings;

  ThrustVectorOCPSettings<Scalar> ocpSettings;
};

template <typename Scalar>
class HoverInitializer final
    : public Initializer<Scalar, STATE_DIM, INPUT_DIM> {
 public:
  explicit HoverInitializer(const Scalar alpha = Scalar(1))
      : alpha_(alpha), oneMinusAlpha_(Scalar(1) - alpha) {
    assert(alpha_ >= Scalar(0));
    assert(alpha_ <= Scalar(1));
  }

  void compute(const Scalar time, const Vector<Scalar, STATE_DIM>& state,
               const Scalar nextTime, Vector<Scalar, INPUT_DIM>& input,
               Vector<Scalar, STATE_DIM>& nextState) override {
    (void)time;
    const Scalar dt = nextTime - time;
    input = hoverInput();

    for (int i = 0; i < 3; ++i) {
      const Scalar commandAcceleration = state(i + 6) + input(i);
      const Scalar effectiveAcceleration =
          oneMinusAlpha_ * state(i + 3) + alpha_ * commandAcceleration;
      nextState(i) = state(i) + dt * effectiveAcceleration;
      nextState(i + 3) = effectiveAcceleration;
      nextState(i + 6) = commandAcceleration;
    }
  }

  static Vector<Scalar, INPUT_DIM> hoverInput() {
    return Vector<Scalar, INPUT_DIM>{Scalar(0.0), Scalar(0.0), Scalar(0.0)};
  }

 private:
  Scalar alpha_;
  Scalar oneMinusAlpha_;
};

template <typename Scalar, size_t PredictLength>
class ThrustVectorOptimalControlProblem {
 public:
  using Problem_t = ThrustVectorProblem<Scalar, PredictLength>;
  using TrackCost_t =
      ThrustVectorTrackCost<Scalar, static_cast<int>(PredictLength + 1)>;
  using DirectionCost_t =
      ThrustDirectionChangeCost<Scalar, static_cast<int>(PredictLength + 1)>;
  using FinalCost_t =
      ThrustVectorTrackFinalCost<Scalar, static_cast<int>(PredictLength + 1)>;
  using ThrustVectorDynamicSystem_t = ThrustVectorDynamicSystem<Scalar>;

  ThrustVectorOptimalControlProblem(
      const ThrustVectorOCPSettings<Scalar>& settings)
      : trackCost_(settings.Q, settings.R, 0),
        directionCost_(settings.alpha, settings.weight, 1),
        finalCost_(settings.Qf, 0),
        dynamics_(settings.alpha) {
    problem_.dynamicsPtr = &dynamics_;
    problem_.cost.add(trackCost_);
    problem_.cost.add(directionCost_);
    problem_.finalCost.add(finalCost_);
  }
  ~ThrustVectorOptimalControlProblem() = default;

  Problem_t& problem() { return problem_; }
  const Problem_t& problem() const { return problem_; }

 protected:
  TrackCost_t trackCost_;
  DirectionCost_t directionCost_;
  FinalCost_t finalCost_;
  ThrustVectorDynamicSystem_t dynamics_;
  Problem_t problem_;
};

template <typename Scalar, size_t PredictLength>
class ThrustVectorILQR
    : public ThrustVectorOptimalControlProblem<Scalar, PredictLength> {
 public:
  using Descriptor_t = iLQRDescriptor<
      Scalar, TranscriptionConfig<Dimensions<STATE_DIM, INPUT_DIM>,
                                  Horizon<PredictLength>, DiscreteDynamics>>;
  using Solver_t = iLQR<Descriptor_t>;
  using HoverInitializer_t = HoverInitializer<Scalar>;
  using StateVector_t = typename Solver_t::StateVector_t;
  using InputVector_t = typename Solver_t::InputVector_t;
  using TimeTrajectory_t = typename Solver_t::TimeTrajectory_t;
  using StateTrajectory_t = typename Solver_t::StateTrajectory_t;
  using InputTrajectory_t = typename Solver_t::InputTrajectory_t;

  explicit ThrustVectorILQR(const ThrustVectorILQRSettings<Scalar>& settings)
      : ThrustVectorOptimalControlProblem<Scalar, PredictLength>(
            settings.ocpSettings),
        initializer_(settings.ocpSettings.alpha),
        solver_(settings.ddpSettings, this->problem_, &initializer_) {}

  ~ThrustVectorILQR() = default;
  ThrustVectorILQR(const ThrustVectorILQR& other) = delete;
  ThrustVectorILQR& operator=(const ThrustVectorILQR& rhs) = delete;
  ThrustVectorILQR(ThrustVectorILQR&& other) = delete;
  ThrustVectorILQR& operator=(ThrustVectorILQR&& rhs) = delete;

  Solver_t& solver() { return solver_; }
  static InputVector_t hoverInput() {
    return InputVector_t{Scalar(0.0), Scalar(0.0), Scalar(0.0)};
  }

  void setDesireTrajectory(const float current_time,
                           const Vector<Scalar, 3>& vel_sp,
                           const Vector<Scalar, 3>& last_input) {
    auto& targetTrajectory = solver_.targetTrajectory();

    for (size_t i = 0; i < PredictLength + 1; ++i) {
      targetTrajectory.timeTrajectory[i] =
          current_time + static_cast<float>(i) * solver_.ddpSettings().timeStep;
      targetTrajectory.stateTrajectory[i].setZero();
      targetTrajectory.stateTrajectory[i].template head<3>() = vel_sp;
      if (i == 0) {
        targetTrajectory.stateTrajectory[i].template segment<3>(6) = last_input;
        targetTrajectory.inputTrajectory[i] = -last_input;
      } else {
        targetTrajectory.inputTrajectory[i] = hoverInput();
      }
    }
  }

 private:
  HoverInitializer_t initializer_;
  Solver_t solver_;
};

template <typename Scalar, size_t PredictLength>
inline ThrustVectorILQR<Scalar, PredictLength> createThrustVectorProblem(
    const ThrustVectorILQRSettings<Scalar>& settings) {
  return ThrustVectorILQR<Scalar, PredictLength>(settings);
}
}  // namespace thrust_vector_lag
