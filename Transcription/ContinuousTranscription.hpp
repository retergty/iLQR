/**
 * @file ContinuousTranscription.hpp
 * @brief 连续时间模型转录：将连续动力学与 running cost 转换为离散 LQ 节点。
 */
#pragma once

#include "LinearQuadraticApproximator.hpp"
#include "ModelData.hpp"
#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "SensitivityIntegrator.hpp"
#include "iLQRDescriptorTraits.hpp"

/**
 * @brief 连续时间转录器。
 *
 * 该类保持当前求解器的语义：动力学通过灵敏度积分器离散化，运行代价按
 * 时间步长缩放后进入离散 Riccati 递推。
 * @tparam Descriptor iLQR 描述类型。
 */
template <typename Descriptor>
class ContinuousTranscription {
 public:
  using Traits = iLQRDescriptorTraits<Descriptor>;
  using Scalar = typename Traits::Scalar;
  using TranscriptionConfig = typename Traits::TranscriptionConfig;
  using ConstraintConfig = typename Traits::ConstraintConfig;

  static constexpr int XDim = Traits::XDim;
  static constexpr int UDim = Traits::UDim;

  using OptimalControlProblem_t =
      OptimalControlProblem<Scalar, TranscriptionConfig, ConstraintConfig>;
  using TargetTrajectories_t = typename Traits::TargetTrajectories_t;
  using StateVector_t = typename Traits::StateVector_t;
  using InputVector_t = typename Traits::InputVector_t;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;
  using Approximator_t = LinearQuadraticApproximator<Descriptor>;
  using DynamicsDiscretizer_t = EK2DynamicsDiscretizer<Scalar, XDim, UDim>;
  using IntermediateMultiplierCollection_t = MultiplierCollection<
      Scalar, typename Traits::IntermediateStageConstraintLayout_t>;
  using FinalMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           typename Traits::FinalStageConstraintLayout_t>;

  ContinuousTranscription() = default;
  ~ContinuousTranscription() = default;

  /**
   * @brief 生成中间节点的离散 LQ 近似。
   *
   * 连续模式下，代价函数近似被解释为 running cost density，因此需要乘以
   * dt；动力学先由连续系统离散化，再转成无仿射缺陷的偏差动力学。
   */
  void approximateIntermediate(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, Scalar time,
      const StateVector_t& state, const InputVector_t& input, Scalar timeStep,
      const IntermediateMultiplierCollection_t& multipliers,
      ModelData_t& modelData) {
    const ModelData_t continuousTimeModelData =
        Approximator_t::approximateIntermediateCostLQ(
            problem, targetTrajectories, time, state, input, multipliers);

    modelData.time = continuousTimeModelData.time;
    modelData.dynamics = discretizer_.sensitivityDiscretize(
        *problem.dynamicsPtr, time, state, input, timeStep);
    modelData.dynamics.f.setZero();
    modelData.cost = continuousTimeModelData.cost * timeStep;
  }

  /**
   * @brief 生成终端节点 LQ 近似。
   *
   * 终端代价不属于 running cost，不做 dt 缩放。
   */
  void approximateFinal(const OptimalControlProblem_t& problem,
                        const TargetTrajectories_t& targetTrajectories,
                        Scalar time, const StateVector_t& state,
                        const FinalMultiplierCollection_t& multipliers,
                        ModelData_t& modelData) const {
    Approximator_t::approximateFinalLQ(problem, targetTrajectories, time, state,
                                       multipliers, modelData);
  }

 private:
  DynamicsDiscretizer_t discretizer_;
};
