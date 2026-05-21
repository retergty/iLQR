#pragma once
#include <array>

#include "ControlledSystemBase.hpp"
#include "DDPData.hpp"
#include "DiscreteTimeRiccatiEquations.hpp"
#include "DualSolution.hpp"
#include "Initializer.hpp"
#include "InitializerRollout.hpp"
#include "LinearController.hpp"
#include "LinearQuadraticApproximator.hpp"
#include "Metrics.hpp"
#include "ModelData.hpp"
#include "OptimalControlProblem.hpp"
#include "PerformanceIndex.hpp"
#include "PrimalSolution.hpp"
#include "ProblemMetrics.hpp"
#include "QuadraticApproximation.hpp"
#include "RiccatiModification.hpp"
#include "RolloutBase.hpp"
#include "SearchStrategyBase.hpp"
#include "SensitivityIntegrator.hpp"
#include "TimeTriggeredRollout.hpp"
#include "Types.hpp"
#include "iLQRDescriptorTraits.hpp"

template <typename Descriptor>
struct iLQRTypes {
  using Traits = iLQRDescriptorTraits<Descriptor>;
  using Descriptor_t = typename Traits::Descriptor_t;
  using Scalar = typename Traits::Scalar;
  using TranscriptionConfig = typename Traits::TranscriptionConfig;
  using Dims = typename Traits::Dims;
  using Horizon = typename Traits::Horizon;
  using ConstraintConfig = typename Traits::ConstraintConfig;

  static constexpr int XDim = Traits::XDim;
  static constexpr int UDim = Traits::UDim;
  static constexpr std::size_t PredictLength = Traits::PredictLength;

  static constexpr int StateEq = Traits::StateEq;
  static constexpr int StateIneq = Traits::StateIneq;
  static constexpr int StateInputEq = Traits::StateInputEq;
  static constexpr int StateInputIneq = Traits::StateInputIneq;
  static constexpr int FinalStateEq = Traits::FinalStateEq;
  static constexpr int FinalStateIneq = Traits::FinalStateIneq;

  using StateVector_t = Vector<Scalar, XDim>;
  using InputVector_t = Vector<Scalar, UDim>;
  using LvVector_t = Vector<Scalar, UDim>;
  using KmMatrix_t = Matrix<Scalar, UDim, XDim>;
  using SmMatrix_t = Matrix<Scalar, XDim, XDim>;
  using SvVector_t = Vector<Scalar, XDim>;
  using GmMatrix_t = Matrix<Scalar, UDim, XDim>;
  using HmMatrix_t = Matrix<Scalar, UDim, UDim>;
  using GvVector_t = Vector<Scalar, UDim>;

  using TimeTrajectory_t = std::array<Scalar, PredictLength + 1>;
  using StateTrajectory_t = std::array<StateVector_t, PredictLength + 1>;
  using InputTrajectory_t = std::array<InputVector_t, PredictLength + 1>;

  using ControlledSystemBase_t = ControlledSystemBase<Scalar, XDim, UDim>;
  using SystemDynamicsBase_t = SystemDynamicsBase<Scalar, XDim, UDim>;
  using RolloutBase_t = RolloutBase<Scalar, XDim, UDim>;
  using InitializerRollout_t = InitializerRollout<Scalar, XDim, UDim>;
  using Initializer_t = Initializer<Scalar, XDim, UDim>;
  using TimeTriggeredRollout_t = TimeTriggeredRollout<Scalar, XDim, UDim>;
  using RolloutTrajectoryPointer_t =
      RolloutTrajectoryPointer<Scalar, XDim, UDim>;

  using OptimalControlProblem_t =
      OptimalControlProblem<Scalar, TranscriptionConfig, ConstraintConfig>;
  using LinearController_t =
      LinearController<Scalar, XDim, UDim, PredictLength + 1>;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;
  using IntermediateMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           typename Traits::IntermediateStageConstraintLayout_t>;
  using FinalMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           typename Traits::FinalStageConstraintLayout_t>;
  using IntermediateMetrics_t =
      Metrics<Scalar, Dims, typename Traits::IntermediateStageConstraintLayout_t>;
  using FinalMetrics_t =
      Metrics<Scalar, Dims, typename Traits::FinalStageConstraintLayout_t>;
  using PrimalSolution_t = PrimalSolution<Scalar, TranscriptionConfig>;
  using DualSolution_t = DualSolution<Scalar, Horizon, ConstraintConfig>;
  using DualSolutionRef_t = DualSolutionRef<Scalar, Horizon, ConstraintConfig>;
  using LinearQuadraticApproximator_t = LinearQuadraticApproximator<Descriptor>;
  using PrimalDataContainer_t =
      PrimalDataContainer<Scalar, TranscriptionConfig, ConstraintConfig>;
  using DualDataContainer_t =
      DualDataContainer<Scalar, TranscriptionConfig, ConstraintConfig>;
  using ProblemMetrics_t =
      ProblemMetrics<Scalar, TranscriptionConfig, ConstraintConfig>;
  using PerformanceIndex_t = PerformanceIndex<Scalar>;
  using IntermediatePerformanceIndexTrajectory_t =
      std::array<PerformanceIndex_t, PredictLength>;
  using IntermediateMultiplierTrajectory_t =
      std::array<IntermediateMultiplierCollection_t, PredictLength>;
  using ModelDataTrajectory_t = std::array<ModelData_t, PredictLength>;
  using SearchStrategySolution_t = SearchStrategySolution<Descriptor>;
  using SearchStrategySolutionRef_t = SearchStrategySolutionRef<Descriptor>;
  using SearchStrategyBase_t = SearchStrategyBase<Descriptor>;
  using EK2DynamicsDiscretizer_t = EK2DynamicsDiscretizer<Scalar, XDim, UDim>;
  using ValueFunctionQuadraticApproximation_t =
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>;
  using ValueFunctionTrajectory_t =
      std::array<ValueFunctionQuadraticApproximation_t, PredictLength + 1>;
  using DiscreteTimeRiccatiEquations_t =
      DiscreteTimeRiccatiEquations<Scalar, XDim, UDim>;
  using RiccatiModification_t = RiccatiModification<Scalar, XDim, UDim>;
};