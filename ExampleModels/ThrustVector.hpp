#pragma once

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
 * 该示例定义了一个用于速度外环的离散最优控制问题。
 *
 * 当前模型不是完整的矢量推力无人机刚体动力学，而是在预测时域内
 * 姿态已知且近似不变时的 NED 平动增量 MPC 模型：
 * - 状态 x = [v, a_prev]，其中 v 为 NED 速度，a_prev 为上一拍实际生效的
 *   NED 总加速度（应写入实际发送或估计值，而非限幅前期望指令）。
 * - 输入 u 为本拍 NED 总加速度增量 Δa_tot，不是机体系推力。
 * - 动力学：a_tot,k = a_prev + u，
 *           v_{k+1} = v_k + dt * a_tot,k，
 *           a_prev_{k+1} = a_tot,k。
 *
 * 推力方向变化代价先从总加速度还原 a_tot，再减去重力得到推力加速度，
 * 惩罚相邻两拍推力方向变化；低推力区通过门控权重减弱该代价。
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
      approximation.setZero();
      dt = 0;
      dirty = true;
      A.template topLeftCorner<3, 3>().setIdentity();
      A.template bottomRightCorner<3, 3>().setIdentity();
      approximation.dfdx = A;
      approximation.dfdu.template bottomRows<3>().setIdentity();
    }
    Matrix<Scalar, STATE_DIM, STATE_DIM> A;
    LinearApproximation_t approximation;
    Scalar dt;
    bool dirty;
  };

  ThrustVectorDynamicSystem() = default;
  ~ThrustVectorDynamicSystem() override = default;

  void updateCache(const Scalar dt) {
    if (cache_.dirty ||
        std::abs(cache_.dt - dt) > std::numeric_limits<Scalar>::epsilon()) {
      for (int i = 0; i < 3; ++i) {
        cache_.A(i, i + 3) = dt;
        cache_.approximation.dfdu(i, i) = dt;
      }
      cache_.approximation.dfdx = cache_.A;
      cache_.dt = dt;
      cache_.dirty = false;
    }
  }

  void linearApproximation(Scalar t, const Vector<Scalar, STATE_DIM>& x,
                           const Vector<Scalar, INPUT_DIM>& u, Scalar dt,
                           LinearApproximation_t& approximation) override {
    approximation.f = computeMap(t, x, u, dt);
    approximation.dfdx = cache_.A;
    approximation.dfdu = cache_.approximation.dfdu;
  }

  void deviationLinearApproximation(
      Scalar t, const Vector<Scalar, STATE_DIM>& x,
      const Vector<Scalar, INPUT_DIM>& u, Scalar dt,
      LinearApproximation_t& approximation) override {
    (void)t;
    (void)x;
    (void)u;

    updateCache(dt);
    approximation = cache_.approximation;
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
    for (int i = 0; i < 3; ++i) {
      next_state(i) = x(i) + dt * currentAcceleration(i);
      next_state(i + 3) = currentAcceleration(i);
    }

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
        R_(R) {}
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
    const Vector<Scalar, INPUT_DIM> inputDeviation =
        input - interpolateInputReference(indexAlpha, inputTrajectory);
    const Vector<Scalar, 3> weightedVelocityDeviation = Qv_ * velocityDeviation;
    const Vector<Scalar, INPUT_DIM> weightedInputDeviation =
        R_ * inputDeviation;

    return Scalar(0.5) * velocityDeviation.dot(weightedVelocityDeviation) +
           Scalar(0.5) * inputDeviation.dot(weightedInputDeviation);
  }

  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory,
      ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>&
          addAppro) override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    const Vector<Scalar, 3> velocityDeviation =
        state.template head<3>() -
        interpolateVelocityReference(indexAlpha, stateTrajectoy);
    const Vector<Scalar, INPUT_DIM> inputDeviation =
        input - interpolateInputReference(indexAlpha, inputTrajectory);
    const Vector<Scalar, 3> weightedVelocityDeviation = Qv_ * velocityDeviation;
    const Vector<Scalar, INPUT_DIM> weightedInputDeviation =
        R_ * inputDeviation;

    addAppro.f +=
        Scalar(0.5) * velocityDeviation.dot(weightedVelocityDeviation) +
        Scalar(0.5) * inputDeviation.dot(weightedInputDeviation);
    addAppro.dfdx.template head<3>() += weightedVelocityDeviation;
    addAppro.dfdu += weightedInputDeviation;
    addAppro.dfdxx.template topLeftCorner<3, 3>() += Qv_;
    addAppro.dfduu += R_;
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
};

template <typename Scalar, int ArrayLength>
class ThrustVectorDiagonalTrackCost final
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {
 public:
  ThrustVectorDiagonalTrackCost(const Matrix<Scalar, STATE_DIM, STATE_DIM>& Q,
                                const Matrix<Scalar, INPUT_DIM, INPUT_DIM>& R,
                                int cost_number)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number) {
    for (int i = 0; i < 3; ++i) {
      QvDiagonal_(i) = Q(i, i);
      RDiagonal_(i) = R(i, i);
    }
  }
  ~ThrustVectorDiagonalTrackCost() override = default;

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    Scalar value = 0;
    for (int i = 0; i < 3; ++i) {
      const Scalar velocityDeviation =
          state(i) -
          interpolateStateReferenceComponent(indexAlpha, stateTrajectoy, i);
      const Scalar inputDeviation =
          input(i) -
          interpolateInputReferenceComponent(indexAlpha, inputTrajectory, i);

      value += QvDiagonal_(i) * velocityDeviation * velocityDeviation +
               RDiagonal_(i) * inputDeviation * inputDeviation;
    }
    return Scalar(0.5) * value;
  }

  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory,
      ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>&
          addAppro) override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    Scalar value = 0;
    for (int i = 0; i < 3; ++i) {
      const Scalar velocityDeviation =
          state(i) -
          interpolateStateReferenceComponent(indexAlpha, stateTrajectoy, i);
      const Scalar inputDeviation =
          input(i) -
          interpolateInputReferenceComponent(indexAlpha, inputTrajectory, i);

      value += QvDiagonal_(i) * velocityDeviation * velocityDeviation +
               RDiagonal_(i) * inputDeviation * inputDeviation;

      addAppro.dfdx(i) += QvDiagonal_(i) * velocityDeviation;
      addAppro.dfdu(i) += RDiagonal_(i) * inputDeviation;
      addAppro.dfdxx(i, i) += QvDiagonal_(i);
      addAppro.dfduu(i, i) += RDiagonal_(i);
    }
    addAppro.f += Scalar(0.5) * value;
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

  Vector<Scalar, 3> QvDiagonal_;
  Vector<Scalar, INPUT_DIM> RDiagonal_;
};

template <typename Scalar, int ArrayLength>
class ThrustVectorTrackFinalCost final
    : public StateCost<Scalar, STATE_DIM, ArrayLength> {
 public:
  ThrustVectorTrackFinalCost(const Matrix<Scalar, STATE_DIM, STATE_DIM>& Q,
                             int cost_number)
      : StateCost<Scalar, STATE_DIM, ArrayLength>(cost_number),
        Qv_(Q.template topLeftCorner<3, 3>()) {}
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

  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, 0>& addAppro)
      override {
    const Vector<Scalar, 3> velocityDeviation =
        state.template head<3>() -
        interpolateVelocityReference(time, timeTrajectory, stateTrajectory);
    const Vector<Scalar, 3> weightedVelocityDeviation = Qv_ * velocityDeviation;

    addAppro.f +=
        Scalar(0.5) * velocityDeviation.dot(weightedVelocityDeviation);
    addAppro.dfdx.template head<3>() += weightedVelocityDeviation;
    addAppro.dfdxx.template topLeftCorner<3, 3>() += Qv_;
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
};

template <typename Scalar, int ArrayLength>
class ThrustVectorDiagonalTrackFinalCost final
    : public StateCost<Scalar, STATE_DIM, ArrayLength> {
 public:
  ThrustVectorDiagonalTrackFinalCost(
      const Matrix<Scalar, STATE_DIM, STATE_DIM>& Q, int cost_number)
      : StateCost<Scalar, STATE_DIM, ArrayLength>(cost_number) {
    for (int i = 0; i < 3; ++i) {
      QvDiagonal_(i) = Q(i, i);
    }
  }
  ~ThrustVectorDiagonalTrackFinalCost() override = default;

  Scalar getValue(Scalar time, const Vector<Scalar, STATE_DIM>& state,
                  const std::array<Scalar, ArrayLength>& timeTrajectory,
                  const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>&
                      stateTrajectory) override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    Scalar value = 0;
    for (int i = 0; i < 3; ++i) {
      const Scalar velocityDeviation =
          state(i) -
          interpolateStateReferenceComponent(indexAlpha, stateTrajectory, i);
      value += QvDiagonal_(i) * velocityDeviation * velocityDeviation;
    }
    return Scalar(0.5) * value;
  }

  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, 0>& addAppro)
      override {
    const auto indexAlpha =
        LinearInterpolation::timeSegment(time, timeTrajectory);
    Scalar value = 0;
    for (int i = 0; i < 3; ++i) {
      const Scalar velocityDeviation =
          state(i) -
          interpolateStateReferenceComponent(indexAlpha, stateTrajectory, i);
      value += QvDiagonal_(i) * velocityDeviation * velocityDeviation;

      addAppro.dfdx(i) += QvDiagonal_(i) * velocityDeviation;
      addAppro.dfdxx(i, i) += QvDiagonal_(i);
    }
    addAppro.f += Scalar(0.5) * value;
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

  Vector<Scalar, 3> QvDiagonal_;
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
      : ThrustDirectionChangeCost(Weight, cost_number) {}

  ThrustDirectionChangeCost(Scalar weight, int cost_number)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number),
        weight_(weight),
        gravityVector_{Scalar(0), Scalar(0), GravityAcceleration} {}

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
        state.template tail<3>() - gravityVector_;
    const Vector<Scalar, 3> thrustAcceleration =
        state.template tail<3>() + input - gravityVector_;
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

  /** @brief 计算并累加代价的二次近似（状态-输入）。 */
  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory,
      ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>&
          addAppro) override {
    (void)time;
    (void)timeTrajectory;
    (void)inputTrajectory;
    (void)stateTrajectory;

    const Vector<Scalar, 3> lastThrustAcceleration =
        state.template tail<3>() - gravityVector_;
    const Vector<Scalar, 3> thrustAcceleration =
        state.template tail<3>() + input - gravityVector_;
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
    const Matrix<Scalar, 3, 3> stateAccelerationJacobian =
        thrustAccelerationJacobian - lastThrustAccelerationJacobian;

    const Scalar gate =
        lowAccelerationGate(lastThrustAcceleration, thrustAcceleration);
    const Scalar effectiveWeight = gate * weight_;

    addAppro.f += Scalar(0.5) * effectiveWeight * rk.dot(rk);
    addAppro.dfdx.template tail<3>() +=
        effectiveWeight * stateAccelerationJacobian.transpose() * rk;
    addAppro.dfdu +=
        effectiveWeight * thrustAccelerationJacobian.transpose() * rk;
    addAppro.dfdxx.template bottomRightCorner<3, 3>() +=
        effectiveWeight * stateAccelerationJacobian.transpose() *
        stateAccelerationJacobian;
    addAppro.dfduu += effectiveWeight * thrustAccelerationJacobian.transpose() *
                      thrustAccelerationJacobian;
    addAppro.dfdux.template topRightCorner<3, 3>() +=
        effectiveWeight * thrustAccelerationJacobian.transpose() *
        stateAccelerationJacobian;
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

  const Scalar weight_;
  const Vector<Scalar, 3> gravityVector_;
};

template <typename Scalar>
struct ThrustVectorOCPSettings {
  // 权重
  Matrix<Scalar, STATE_DIM, STATE_DIM> Q;
  Matrix<Scalar, INPUT_DIM, INPUT_DIM> R;
  Matrix<Scalar, STATE_DIM, STATE_DIM> Qf;
  Scalar weight;
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
  HoverInitializer() {};
  void compute(const Scalar time, const Vector<Scalar, STATE_DIM>& state,
               const Scalar nextTime, Vector<Scalar, INPUT_DIM>& input,
               Vector<Scalar, STATE_DIM>& nextState) override {
    (void)time;
    const Scalar dt = nextTime - time;
    input = hoverInput();
    const Vector<Scalar, INPUT_DIM> currentAcceleration =
        state.template tail<INPUT_DIM>() + input;

    nextState.template head<3>() =
        state.template head<3>() + dt * currentAcceleration;
    nextState.template tail<3>() = currentAcceleration;
  }
  static Vector<Scalar, INPUT_DIM> hoverInput() {
    return Vector<Scalar, INPUT_DIM>{Scalar(0.0), Scalar(0.0), Scalar(0.0)};
  }
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
        directionCost_(settings.weight, 1),
        finalCost_(settings.Qf, 0) {
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
        targetTrajectory.stateTrajectory[i].template tail<3>() = last_input;
        targetTrajectory.inputTrajectory[i] = -last_input;
      } else {
        targetTrajectory.stateTrajectory[i].template tail<3>().setZero();
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
}  // namespace thrust_vector
