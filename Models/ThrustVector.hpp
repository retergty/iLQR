#pragma once

#include <limits>

#include "Cost.hpp"
#include "DiscreteSystemDynamicsBase.hpp"
#include "LinearApproximation.hpp"
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
      A(0, 0) = 1;
      A(1, 1) = 1;
      A(2, 2) = 1;
      B(3, 0) = 1;
      B(4, 1) = 1;
      B(5, 2) = 1;
      approximation.dfdx = A;
    }
    Matrix<Scalar, STATE_DIM, STATE_DIM> A;
    Matrix<Scalar, STATE_DIM, INPUT_DIM> B;
    Vector<Scalar, STATE_DIM> bias;
    LinearApproximation_t approximation;
    Scalar dt;
  };

  ThrustVectorDynamicSystem(const Scalar mass) : mass_(mass) {
    rotB2w_.setIdentity();
  }
  ~ThrustVectorDynamicSystem() override = default;

  void updateRotationMatrix(const Matrix<Scalar, 3, 3>& rotation) {
    rotB2w_ = rotation;
    cache_.dt = 0;
  }

  void updateCache(const Scalar dt) {
    if (std::abs(cache_.dt - dt) > std::numeric_limits<Scalar>::epsilon()) {
      cache_.B(0, 0) = rotB2w_(0, 0) / mass_ * dt;
      cache_.B(0, 1) = rotB2w_(0, 1) / mass_ * dt;
      cache_.B(0, 2) = rotB2w_(0, 2) / mass_ * dt;
      cache_.B(1, 0) = rotB2w_(1, 0) / mass_ * dt;
      cache_.B(1, 1) = rotB2w_(1, 1) / mass_ * dt;
      cache_.B(1, 2) = rotB2w_(1, 2) / mass_ * dt;
      cache_.B(2, 0) = rotB2w_(2, 0) / mass_ * dt;
      cache_.B(2, 1) = rotB2w_(2, 1) / mass_ * dt;
      cache_.B(2, 2) = rotB2w_(2, 2) / mass_ * dt;

      cache_.bias(2, 0) = -dt * Gravity_Acc;
      cache_.dt = dt;
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
    (void)x;
    (void)u;

    Vector<Scalar, STATE_DIM> next_state;
    next_state(0) = x(0);
    next_state(1) = x(1);
    next_state(2) = x(2);
    next_state(3) = 0;
    next_state(4) = 0;
    next_state(5) = 0;

    updateCache(dt);

    next_state += cache_.B * u + cache_.bias;

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
                        const Matrix<Scalar, INPUT_DIM, INPUT_DIM>& R)
      : QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength>(Q,
                                                                           R) {}
  ~ThrustVectorTrackCost() override = default;

 private:
  ThrustVectorTrackCost(const ThrustVectorTrackCost& other) = default;
};

template <typename Scalar, int ArrayLength>
class ThrustDirectionChangeCost final
    : public StateInputCost<Scalar, STATE_DIM, INPUT_DIM, ArrayLength> {};
};  // namespace thrust_vector