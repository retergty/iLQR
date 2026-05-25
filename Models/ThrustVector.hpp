#pragma once

#include <cmath>
#include <limits>

#include "Cost.hpp"
#include "DiscreteSystemDynamicsBase.hpp"
#include "OptimalControlProblem.hpp"
#include "QuadraticStateCost.hpp"
#include "Types.hpp"
#include "iLQRDescriptor.hpp"
/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
/**
 * 该示例定义了一个最优控制问题，其中运动学建模的
 * 是一个矢量推力无人机，状态空间是6维
 * 全局坐标系三轴速度+上一拍机体系推力，输入是3维机体系推力
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
  static constexpr Scalar Gravity_Acc = 9.8;
  using LinearApproximation_t =
      typename DiscreteSystemDynamicsBase<Scalar, STATE_DIM,
                                          INPUT_DIM>::LinearApproximation_t;

  struct PreCompCache {
    PreCompCache() {
      A.setZero();
      B.setZero();
      bias.setZero();
      approximation.setZero();
      dt = 0;
      dirty = true;
      A.template topLeftCorner<3, 3>().setIdentity();
      B.template bottomRows<3>().setIdentity();
      approximation.dfdx = A;
    }
    Matrix<Scalar, STATE_DIM, STATE_DIM> A;
    Matrix<Scalar, STATE_DIM, INPUT_DIM> B;
    Vector<Scalar, STATE_DIM> bias;
    LinearApproximation_t approximation;
    Scalar dt;
    bool dirty;
  };

  ThrustVectorDynamicSystem(const Scalar mass) : mass_(mass) {
    rotB2w_.setIdentity();
  }
  ~ThrustVectorDynamicSystem() override = default;

  void updateRotationMatrix(const Matrix<Scalar, 3, 3>& rotation) {
    rotB2w_ = rotation;
    cache_.dirty = true;
  }

  void updateCache(const Scalar dt) {
    if (cache_.dirty ||
        std::abs(cache_.dt - dt) > std::numeric_limits<Scalar>::epsilon()) {
      cache_.B.template topRows<3>() = (dt / mass_) * rotB2w_;
      cache_.bias(2) = -dt * Gravity_Acc;
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
    next_state.template head<3>() = x.template head<3>() +
                                    cache_.B.template topRows<3>() * u +
                                    cache_.bias.template head<3>();
    next_state.template tail<3>() = u;

    return next_state;
  }

 private:
  Scalar mass_;
  Matrix<Scalar, 3, 3> rotB2w_;
  PreCompCache cache_;
};

template <typename Scalar, int ArrayLength>
class ThrustVectorTrackCost final
    : public QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM,
                                     ArrayLength> {
 public:
  ThrustVectorTrackCost(const Matrix<Scalar, STATE_DIM, STATE_DIM>& Q,
                        const Matrix<Scalar, INPUT_DIM, INPUT_DIM>& R,
                        int cost_number)
      : QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(
            Q, R, cost_number) {}
  ~ThrustVectorTrackCost() override = default;

 private:
  ThrustVectorTrackCost(const ThrustVectorTrackCost& other) = default;
};

template <typename Scalar, int ArrayLength>
class ThrustDirectionChangeCost final
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {
  /** @brief 获取代价值。 */
 public:
  static constexpr Scalar epsilon = 1e-4;
  static constexpr Scalar Weight = 1;

  explicit ThrustDirectionChangeCost(int cost_number = 0)
      : StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(cost_number) {
  }

  Scalar getValue(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      const override {
    (void)time;
    (void)timeTrajectory;
    (void)inputTrajectory;
    (void)stateTrajectoy;

    const Vector<Scalar, 3> last_thr{state(3), state(4), state(5)};
    const Vector<Scalar, 3> thr = input;
    Vector<Scalar, 3> last_thr_dir =
        last_thr / std::sqrt(last_thr.dot(last_thr) + epsilon);
    Vector<Scalar, 3> thr_dir = thr / std::sqrt(thr.dot(thr) + epsilon);
    const Vector<Scalar, 3> rk = thr_dir - last_thr_dir;
    return Scalar(0.5) * Weight * rk.dot(rk);
  }

  /** @brief 获取代价的二次近似（状态-输入）。 */
  ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, STATE_DIM>& state,
      const Vector<Scalar, INPUT_DIM>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, STATE_DIM>, ArrayLength>& stateTrajectory,
      const std::array<Vector<Scalar, INPUT_DIM>, ArrayLength>& inputTrajectory)
      const override {
    (void)time;
    (void)timeTrajectory;
    (void)inputTrajectory;
    (void)stateTrajectory;

    ScalarFunctionQuadraticApproximation<Scalar, STATE_DIM, INPUT_DIM> ret;
    ret.setZero();

    const Vector<Scalar, 3> last_thr{state(3), state(4), state(5)};
    const Vector<Scalar, 3> thr = input;
    const Scalar last_thr_norm = std::sqrt(last_thr.dot(last_thr) + epsilon);
    const Scalar thr_norm = std::sqrt(thr.dot(thr) + epsilon);
    const Vector<Scalar, 3> last_thr_dir = last_thr / last_thr_norm;
    const Vector<Scalar, 3> thr_dir = thr / thr_norm;
    const Vector<Scalar, 3> rk = thr_dir - last_thr_dir;

    const Matrix<Scalar, 3, 3> identity = Matrix<Scalar, 3, 3>::Identity();
    const Matrix<Scalar, 3, 3> last_thr_jacobian =
        identity / last_thr_norm -
        (last_thr * last_thr.transpose()) /
            (last_thr_norm * last_thr_norm * last_thr_norm);
    const Matrix<Scalar, 3, 3> thr_jacobian =
        identity / thr_norm -
        (thr * thr.transpose()) / (thr_norm * thr_norm * thr_norm);

    ret.f = Scalar(0.5) * Weight * rk.dot(rk);
    ret.dfdx.template segment<3>(3) =
        -Weight * last_thr_jacobian.transpose() * rk;
    ret.dfdu = Weight * thr_jacobian.transpose() * rk;
    ret.dfdxx.template block<3, 3>(3, 3) =
        Weight * last_thr_jacobian.transpose() * last_thr_jacobian;
    ret.dfduu = Weight * thr_jacobian.transpose() * thr_jacobian;
    ret.dfdux.template block<3, 3>(0, 3) =
        -Weight * thr_jacobian.transpose() * last_thr_jacobian;

    return ret;
  }
};
}  // namespace thrust_vector