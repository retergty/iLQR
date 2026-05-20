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
#include "LineSearchStrategy.hpp"
#include "Types.hpp"

template <typename Descriptor>
struct iLQRTypes {
  using Descriptor_t = Descriptor;
  using Scalar = typename Descriptor::Scalar;
  using Dims = typename Descriptor::Dims;
  using Constraints = typename Descriptor::Constraints;

  static constexpr int XDim = Dims::XDim;
  static constexpr int UDim = Dims::UDim;
  static constexpr std::size_t PredictLength = Dims::PredictLength;

  static constexpr int StateEq = Constraints::StateEq;
  static constexpr int StateIneq = Constraints::StateIneq;
  static constexpr int StateInputEq = Constraints::StateInputEq;
  static constexpr int StateInputIneq = Constraints::StateInputIneq;
  static constexpr int FinalStateEq = Constraints::FinalStateEq;
  static constexpr int FinalStateIneq = Constraints::FinalStateIneq;

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
      OptimalControlProblem<Scalar, XDim, UDim, PredictLength, StateEq,
                            StateInputEq, StateIneq, StateInputIneq,
                            FinalStateEq, FinalStateIneq>;
  using LinearController_t =
      LinearController<Scalar, XDim, UDim, PredictLength + 1>;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;
  using IntermediateMultiplierCollection_t =
      MultiplierCollection<Scalar, StateEq, StateIneq, StateInputEq,
                           StateInputIneq>;
  using FinalMultiplierCollection_t =
      MultiplierCollection<Scalar, FinalStateEq, FinalStateIneq, 0, 0>;
  using IntermediateMetrics_t = Metrics<Scalar, XDim, UDim, StateEq, StateIneq,
                                        StateInputEq, StateInputIneq>;
  using FinalMetrics_t =
      Metrics<Scalar, XDim, UDim, FinalStateEq, FinalStateIneq, 0, 0>;
  using PrimalSolution_t = PrimalSolution<Scalar, XDim, UDim, PredictLength>;
  using DualSolution_t =
      DualSolution<Scalar, StateEq, StateIneq, StateInputEq, StateInputIneq,
                   FinalStateEq, FinalStateIneq, PredictLength>;
  using DualSolutionRef_t =
      DualSolutionRef<Scalar, StateEq, StateIneq, StateInputEq, StateInputIneq,
                      FinalStateEq, FinalStateIneq, PredictLength>;
  using LinearQuadraticApproximator_t =
      LinearQuadraticApproximator<Scalar, XDim, UDim, PredictLength, StateEq,
                                  StateIneq, StateInputEq, StateInputIneq,
                                  FinalStateEq, FinalStateIneq>;
  using PrimalDataContainer_t =
      PrimalDataContainer<Scalar, XDim, UDim, PredictLength, StateEq, StateIneq,
                          StateInputEq, StateInputIneq, FinalStateEq,
                          FinalStateIneq>;
  using DualDataContainer_t =
      DualDataContainer<Scalar, XDim, UDim, PredictLength, StateEq, StateIneq,
                        StateInputEq, StateInputIneq, FinalStateEq,
                        FinalStateIneq>;
  using ProblemMetrics_t =
      ProblemMetrics<Scalar, XDim, UDim, PredictLength, StateEq, StateIneq,
                     StateInputEq, StateInputIneq, FinalStateEq,
                     FinalStateIneq>;
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