/**
 * @file DiscreteTranscription.hpp
 * @brief 离散时间模型转录：直接生成离散 LQ 节点。
 */
#pragma once

#include "LinearQuadraticApproximator.hpp"
#include "ModelData.hpp"
#include "Multiplier.hpp"
#include "OptimalControlProblem.hpp"
#include "iLQRDescriptorTraits.hpp"

/**
 * @brief 离散时间转录器。
 *
 * 该类假设动力学接口已经给出离散状态转移。中间代价默认解释为离散
 * stage cost，因此不再乘以时间步长。
 * @tparam Descriptor iLQR 描述类型。
 */
template <typename Descriptor>
class DiscreteTranscription {
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
  using IntermediateMultiplierCollection_t = MultiplierCollection<
      Scalar, typename Traits::IntermediateStageConstraintLayout_t>;
  using FinalMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           typename Traits::FinalStageConstraintLayout_t>;

  DiscreteTranscription() = default;
  ~DiscreteTranscription() = default;

  /**
   * @brief 生成中间节点的离散 LQ 近似。
   *
   * 离散模式下，动力学直接使用偏差线性化：
   * \f$ \delta x_{k+1}=A_k\delta x_k+B_k\delta u_k \f$。
   * 代价函数近似保持 stage cost 语义，不做 dt 缩放。
   */
  void approximateIntermediate(
      const OptimalControlProblem_t& problem,
      const TargetTrajectories_t& targetTrajectories, Scalar time,
      const StateVector_t& state, const InputVector_t& input, Scalar timeStep,
      const IntermediateMultiplierCollection_t& multipliers,
      ModelData_t& modelData) const {
    const ModelData_t discreteModelData =
        Approximator_t::approximateIntermediateCostLQ(
            problem, targetTrajectories, time, state, input, multipliers);

    modelData.time = discreteModelData.time;
    modelData.dynamics = problem.dynamicsPtr->deviationLinearApproximation(
        time, state, input, timeStep);
    modelData.cost = discreteModelData.cost;
  }

  /**
   * @brief 生成终端节点 LQ 近似。
   */
  void approximateFinal(const OptimalControlProblem_t& problem,
                        const TargetTrajectories_t& targetTrajectories,
                        Scalar time, const StateVector_t& state,
                        const FinalMultiplierCollection_t& multipliers,
                        ModelData_t& modelData) const {
    Approximator_t::approximateFinalLQ(problem, targetTrajectories, time, state,
                                       multipliers, modelData);
  }
};
