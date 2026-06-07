#pragma once

#include <cmath>
#include <limits>

#include "Cost/Cost.hpp"
#include "Dynamics/DiscreteSystemDynamicsBase.hpp"
#include "Misc/LinearInterpolation.hpp"
#include "OptimalControl/OptimalControlProblem.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"
#include "iLQR/iLQRDescriptor.hpp"
/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/**
 * 该示例定义了一个最优控制问题，其中运动学建模的
 * 是一个矢量推力无人机，状态空间是6维
 * NED 全局坐标系三轴速度+上一拍 NED 总加速度，
 * 输入是3维 NED 总加速度增量
 * 同时添加余弦相似度的代价函数
 */

namespace thrust_vector {
static constexpr int STATE_DIM = 6;
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
      A.setZero();
      B.setZero();
      approximation.setZero();
      dt = 0;
      dirty = true;
      A.template topLeftCorner<3, 3>().setIdentity();
      A.template bottomRightCorner<3, 3>().setIdentity();
      B.template bottomRows<3>().setIdentity();
      approximation.dfdx = A;
    }
    Matrix<Scalar, STATE_DIM, STATE_DIM> A;
    Matrix<Scalar, STATE_DIM, INPUT_DIM> B;
    LinearApproximation_t approximation;
    Scalar dt;
    bool dirty;
  };

  ThrustVectorDynamicSystem() = default;
  ~ThrustVectorDynamicSystem() override = default;

  void updateCache(const Scalar dt) {
    if (cache_.dirty ||
        std::abs(cache_.dt - dt) > std::numeric_limits<Scalar>::epsilon()) {
      cache_.A.template topRightCorner<3, 3>() =
          dt * Matrix<Scalar, 3, 3>::Identity();
      cache_.B.template topRows<3>() = dt * Matrix<Scalar, 3, 3>::Identity();
      cache_.approximation.dfdx = cache_.A;
      cache_.dt = dt;
      cache_.dirty = false;
    }
  }

  LinearApproximation_t linearApproximation(Scalar t,
                                            const Vector<Scalar, STATE_DIM>& x,
                                            const Vector<Scalar, INPUT_DIM>& u,
                                            Scalar dt) override {
    LinearApproximation_t approximation;
    approximation.f = computeMap(t, x, u, dt);
    approximation.dfdx = cache_.A;
    approximation.dfdu = cache_.B;
    return approximation;
  }

  LinearApproximation_t deviationLinearApproximation(
      Scalar t, const Vector<Scalar, STATE_DIM>& x,
      const Vector<Scalar, INPUT_DIM>& u, Scalar dt) override {
    (void)t;
    (void)x;
    (void)u;

    updateCache(dt);
    cache_.approximation.dfdu = cache_.B;
    return cache_.approximation;
  }
  Vector<Scalar, STATE_DIM> computeMap(Scalar t,
                                       const Vector<Scalar, STATE_DIM>& x,
                                       const Vector<Scalar, INPUT_DIM>& u,
                                       Scalar dt) override {
    (void)t;

    Vector<Scalar, STATE_DIM> next_state;
    updateCache(dt);
    const Vector<Scalar, INPUT_DIM> currentAcceleration =
        x.template tail<INPUT_DIM>() + u;
    next_state.template head<3>() =
        x.template head<3>() +
        cache_.B.template topRows<3>() * currentAcceleration;
    next_state.template tail<3>() = currentAcceleration;

    return next_state;
  }

 private:
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
    approximation_.dfdxx.template bottomRightCorner<3, 3>() = R_;
    approximation_.dfduu = R_;
    approximation_.dfdux.template topRightCorner<3, 3>() = R_;
  }
  ~ThrustVectorTrackCost() override = default;

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    const Vector<Scalar, 3> velocityDeviation =
        state.template head<3>() -
        interpolateVelocityReference(indexAlpha, stateTrajectoy);
    const Vector<Scalar, INPUT_DIM> accelerationDeviation =
        state.template tail<INPUT_DIM>() + input -
        interpolateInputReference(indexAlpha, inputTrajectory);
    const Vector<Scalar, 3> weightedVelocityDeviation = Qv_ * velocityDeviation;
    const Vector<Scalar, INPUT_DIM> weightedAccelerationDeviation =
        R_ * accelerationDeviation;

    return Scalar(0.5) * velocityDeviation.dot(weightedVelocityDeviation) +
           Scalar(0.5) *
               accelerationDeviation.dot(weightedAccelerationDeviation);
  }

  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    const Vector<Scalar, 3> velocityDeviation =
        state.template head<3>() -
        interpolateVelocityReference(indexAlpha, stateTrajectoy);
    const Vector<Scalar, INPUT_DIM> accelerationDeviation =
        state.template tail<INPUT_DIM>() + input -
        interpolateInputReference(indexAlpha, inputTrajectory);
    const Vector<Scalar, 3> weightedVelocityDeviation = Qv_ * velocityDeviation;
    const Vector<Scalar, INPUT_DIM> weightedAccelerationDeviation =
        R_ * accelerationDeviation;

    approximation_.f =
        Scalar(0.5) * velocityDeviation.dot(weightedVelocityDeviation) +
        Scalar(0.5) * accelerationDeviation.dot(weightedAccelerationDeviation);
    approximation_.dfdx.template head<3>() = weightedVelocityDeviation;
    approximation_.dfdx.template tail<INPUT_DIM>() =
        weightedAccelerationDeviation;
    approximation_.dfdu = weightedAccelerationDeviation;
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
class ThrustInputRateCost final
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {
 public:
  ThrustInputRateCost(const Matrix<Scalar, INPUT_DIM, INPUT_DIM>& S,
                      int cost_number)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number),
        S_(S) {
    approximation_.setZero();
    approximation_.dfduu = S_;
  }

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    (void)time;
    (void)timeTrajectory;
    (void)state;
    (void)stateTrajectoy;
    (void)inputTrajectory;

    return Scalar(0.5) * input.dot(S_ * input);
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
    (void)state;
    (void)stateTrajectory;
    (void)inputTrajectory;

    const Vector<Scalar, INPUT_DIM> weightedInputRate = S_ * input;
    approximation_.f = Scalar(0.5) * input.dot(weightedInputRate);
    approximation_.dfdu = weightedInputRate;
    return approximation_;
  }

 private:
  Matrix<Scalar, INPUT_DIM, INPUT_DIM> S_;
  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
      approximation_;
};

template <typename Scalar, int ArrayLength>
class ThrustDirectionChangeCost final
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {
  /** @brief 获取代价值。 */
 public:
  static constexpr Scalar epsilon = Scalar(1e-4);
  static constexpr Scalar Weight = 1;
  static constexpr Scalar MinThrustAccelerationForDirection = Scalar(1e-2);
  static constexpr Scalar GravityAcceleration = Scalar(9.8);

  explicit ThrustDirectionChangeCost(int cost_number = 0)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number) {
    approximation_.setZero();
  }

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    (void)time;
    (void)timeTrajectory;
    (void)inputTrajectory;
    (void)stateTrajectoy;

    const Vector<Scalar, 3> lastThrustAcceleration =
        state.template tail<3>() - gravityVector();
    const Vector<Scalar, 3> thrustAcceleration =
        state.template tail<3>() + input - gravityVector();
    const Vector<Scalar, 3> lastThrustDirection =
        lastThrustAcceleration /
        std::sqrt(lastThrustAcceleration.dot(lastThrustAcceleration) + epsilon);
    const Vector<Scalar, 3> thrustDirection =
        thrustAcceleration /
        std::sqrt(thrustAcceleration.dot(thrustAcceleration) + epsilon);
    const Vector<Scalar, 3> rk = thrustDirection - lastThrustDirection;
    const Scalar gate =
        lowAccelerationGate(lastThrustAcceleration, thrustAcceleration);
    return Scalar(0.5) * gate * Weight * rk.dot(rk);
  }

  /** @brief 获取代价的二次近似（状态-输入）。 */
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
        state.template tail<3>() - gravityVector();
    const Vector<Scalar, 3> thrustAcceleration =
        state.template tail<3>() + input - gravityVector();
    const Scalar lastThrustAccelerationNorm =
        std::sqrt(lastThrustAcceleration.dot(lastThrustAcceleration) + epsilon);
    const Scalar thrustAccelerationNorm =
        std::sqrt(thrustAcceleration.dot(thrustAcceleration) + epsilon);
    const Vector<Scalar, 3> lastThrustDirection =
        lastThrustAcceleration / lastThrustAccelerationNorm;
    const Vector<Scalar, 3> thrustDirection =
        thrustAcceleration / thrustAccelerationNorm;
    const Vector<Scalar, 3> rk = thrustDirection - lastThrustDirection;

    const Matrix<Scalar, 3, 3> identity = Matrix<Scalar, 3, 3>::Identity();
    const Matrix<Scalar, 3, 3> lastThrustAccelerationJacobian =
        identity / lastThrustAccelerationNorm -
        (lastThrustAcceleration * lastThrustAcceleration.transpose()) /
            (lastThrustAccelerationNorm * lastThrustAccelerationNorm *
             lastThrustAccelerationNorm);
    const Matrix<Scalar, 3, 3> thrustAccelerationJacobian =
        identity / thrustAccelerationNorm -
        (thrustAcceleration * thrustAcceleration.transpose()) /
            (thrustAccelerationNorm * thrustAccelerationNorm *
             thrustAccelerationNorm);
    const Matrix<Scalar, 3, 3> stateAccelerationJacobian =
        thrustAccelerationJacobian - lastThrustAccelerationJacobian;

    const Scalar gate =
        lowAccelerationGate(lastThrustAcceleration, thrustAcceleration);
    const Scalar effectiveWeight = gate * Weight;

    approximation_.f = Scalar(0.5) * effectiveWeight * rk.dot(rk);
    approximation_.dfdx.template tail<3>() =
        effectiveWeight * stateAccelerationJacobian.transpose() * rk;
    approximation_.dfdu =
        effectiveWeight * thrustAccelerationJacobian.transpose() * rk;
    approximation_.dfdxx.template bottomRightCorner<3, 3>() =
        effectiveWeight * stateAccelerationJacobian.transpose() *
        stateAccelerationJacobian;
    approximation_.dfduu = effectiveWeight *
                           thrustAccelerationJacobian.transpose() *
                           thrustAccelerationJacobian;
    approximation_.dfdux.template topRightCorner<3, 3>() =
        effectiveWeight * thrustAccelerationJacobian.transpose() *
        stateAccelerationJacobian;

    return approximation_;
  }

 private:
  static Vector<Scalar, 3> gravityVector() {
    return Vector<Scalar, 3>{Scalar(0), Scalar(0), GravityAcceleration};
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

  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
      approximation_;
};

template <typename Scalar>
Matrix<Scalar, STATE_DIM, STATE_DIM> makeDefaultThrustVectorStateWeight() {
  Matrix<Scalar, STATE_DIM, STATE_DIM> Q;
  Q.setZero();
  Q.template topLeftCorner<3, 3>().setIdentity();
  return Q;
}

template <typename Scalar>
Matrix<Scalar, STATE_DIM, STATE_DIM> makeDefaultThrustVectorFinalWeight() {
  Matrix<Scalar, STATE_DIM, STATE_DIM> Qf;
  Qf.setZero();
  Qf.template topLeftCorner<3, 3>() =
      Scalar(10.0) * Matrix<Scalar, 3, 3>::Identity();
  return Qf;
}

template <typename Scalar>
Matrix<Scalar, INPUT_DIM, INPUT_DIM> makeDefaultThrustVectorInputWeight() {
  return Scalar(1e-3) * Matrix<Scalar, INPUT_DIM, INPUT_DIM>::Identity();
}

template <typename Scalar>
Matrix<Scalar, INPUT_DIM, INPUT_DIM> makeDefaultThrustInputRateWeight() {
  return Scalar(1e-2) * Matrix<Scalar, INPUT_DIM, INPUT_DIM>::Identity();
}

template <typename Scalar, size_t PredictLength>
inline ThrustVectorProblem<Scalar, PredictLength>& createThrustVectorProblem() {
  using Problem_t = ThrustVectorProblem<Scalar, PredictLength>;
  using TrackCost_t =
      ThrustVectorTrackCost<Scalar, static_cast<int>(PredictLength + 1)>;
  using InputRateCost_t =
      ThrustInputRateCost<Scalar, static_cast<int>(PredictLength + 1)>;
  using DirectionCost_t =
      ThrustDirectionChangeCost<Scalar, static_cast<int>(PredictLength + 1)>;
  using FinalCost_t =
      ThrustVectorTrackFinalCost<Scalar, static_cast<int>(PredictLength + 1)>;

  static ThrustVectorDynamicSystem<Scalar> dynamics;
  static const Matrix<Scalar, STATE_DIM, STATE_DIM> Q =
      makeDefaultThrustVectorStateWeight<Scalar>();
  static const Matrix<Scalar, INPUT_DIM, INPUT_DIM> R =
      makeDefaultThrustVectorInputWeight<Scalar>();
  static const Matrix<Scalar, INPUT_DIM, INPUT_DIM> S =
      makeDefaultThrustInputRateWeight<Scalar>();
  static const Matrix<Scalar, STATE_DIM, STATE_DIM> Qf =
      makeDefaultThrustVectorFinalWeight<Scalar>();

  static TrackCost_t trackCost(Q, R, 0);
  static InputRateCost_t inputRateCost(S, 1);
  static DirectionCost_t directionCost(2);
  static FinalCost_t finalCost(Qf, 0);
  static Problem_t problem;
  static bool initialized = false;

  if (!initialized) {
    problem.dynamicsPtr = &dynamics;
    problem.cost.add(trackCost);
    problem.cost.add(inputRateCost);
    problem.cost.add(directionCost);
    problem.finalCost.add(finalCost);
    initialized = true;
  }

  return problem;
}
}  // namespace thrust_vector