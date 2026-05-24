
cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME "tesetCPP")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
project(${PROJECT_NAME} LANGUAGES CXX)

find_package(Eigen3 3.3 CONFIG QUIET)
if(NOT Eigen3_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    eigen
    URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
  )
  FetchContent_MakeAvailable(eigen)

  if(NOT TARGET Eigen3::Eigen AND TARGET eigen)
    add_library(Eigen3::Eigen ALIAS eigen)
  endif()
endif()

add_executable(iLQR iLQR.cpp
  AutomaticDifferentation/FiniteDifferenceMethods.hpp)

set(ILQR_HEADERS
  DDPData.hpp
  DDPSetting.hpp
  HessianCorrection.hpp
  iLQR.hpp
  Types.hpp
  Approximation/ChangeOfInputVariables.hpp
  Approximation/LinearApproximation.hpp
  Approximation/LinearQuadraticApproximator.hpp
  Approximation/QuadraticApproximation.hpp
  AugmentedLagrangian/StateAugmentedLagrangian.hpp
  AugmentedLagrangian/StateAugmentedLagrangianCollection.hpp
  AugmentedLagrangian/StateAugmentedLagrangianInterface.hpp
  AugmentedLagrangian/StateInputAugmentedLagrangian.hpp
  AugmentedLagrangian/StateInputAugmentedLagrangianCollection.hpp
  AugmentedLagrangian/StateInputAugmentedLagrangianInterface.hpp
  AutomaticDifferentation/FiniteDifferenceMethods.hpp
  Constraint/ConstraintOrder.hpp
  Constraint/StateConstraint.hpp
  Constraint/StateInputConstraint.hpp
  Constraint/LinearStateInputConstraint.hpp
  Constraint/LinearStateConstraint.hpp
  Controller/Controller.hpp
  Controller/LinearController.hpp
  Cost/Cost.hpp
  Cost/CostCollection.hpp
  Cost/QuadraticStateCost.hpp
  Dynamics/ControlledSystemBase.hpp
  Dynamics/DiscreteLinearSystemDynamics.hpp
  Dynamics/DiscreteSystemBase.hpp
  Dynamics/DiscreteSystemDynamicsBase.hpp
  Dynamics/DynamicsModeTraits.hpp
  Dynamics/LinearSystemDynamics.hpp
  Dynamics/SystemDynamicsBase.hpp
  Dynamics/SystemDynamicsLinearizer.hpp
  Initialization/DefaultInitializer.hpp
  Initialization/Initializer.hpp
  Initialization/OperatingPoints.hpp
  Integration/Integration.hpp
  Integration/Observer.hpp
  Integration/RungeKuttaDormandPrince5.hpp
  Integration/SensitivityIntegrator.hpp
  Integration/TrapezoidalIntegration.hpp
  IntrusiveList/IntrusiveList.hpp
  Misc/LinearAlgebra.hpp
  Misc/LinearInterpolation.hpp
  Misc/NumericTraits.hpp
  Misc/Numerics.hpp
  OptimalControl/DualSolution.hpp
  OptimalControl/LagrangianMetrics.hpp
  OptimalControl/Metrics.hpp
  OptimalControl/ModelData.hpp
  OptimalControl/Multiplier.hpp
  OptimalControl/OptimalControlProblem.hpp
  OptimalControl/OptimalControlProblemHelperFunction.hpp
  OptimalControl/PerformanceIndex.hpp
  OptimalControl/PrimalSolution.hpp
  OptimalControl/ProblemMetrics.hpp
  OptimalControl/Reference.hpp
  Penalties/AugmentedPenaltyBase.hpp
  Penalties/ModifiedRelaxedBarrierPenalty.hpp
  Penalties/Penalty.hpp
  Penalties/QuadraticPenalty.hpp
  Penalties/SlacknessSquaredHingePenalty.hpp
  Penalties/SmoothAbsolutePenalty.hpp
  Penalties/MultidimensionalPenalty.hpp
  RiccatiEquations/DiscreteTimeRiccatiEquations.hpp
  RiccatiEquations/RiccatiModification.hpp
  Rollout/DiscreteTimeRollout.hpp
  Rollout/InitializerRollout.hpp
  Rollout/RolloutBase.hpp
  Rollout/RolloutTraits.hpp
  Rollout/TimeTriggeredRollout.hpp
  SearchStrategy/LineSearchStrategy.hpp
  SearchStrategy/SearchStrategyBase.hpp
  SearchStrategy/SearchStrategySettings.hpp
  Transcription/ContinuousTranscription.hpp
  Transcription/DiscreteTranscription.hpp
  Transcription/TranscriptionTraits.hpp
  TestTools/QpTrajectories.hpp
  TestTools/QpSolverTypes.hpp
  TestTools/QpSolver.hpp
  TestTools/QpDiscreteTranscription.hpp
  TestTools/Ocs2QpSolver.hpp
  Models/CircularKinematics.hpp
  Models/DoubleIntegratorReachingTask.hpp
  Models/EXP0.hpp
  iLQRDescriptor.hpp
  iLQRDescriptorTraits.hpp
  iLQRTypes.hpp
)

add_library(iLQR_Include INTERFACE)
target_sources(iLQR_Include INTERFACE ${ILQR_HEADERS})
target_include_directories(iLQR_Include INTERFACE Integration Cost Dynamics Controller Approximation OptimalControl Rollout IntrusiveList AugmentedLagrangian Constraint Penalties RiccatiEquations SearchStrategy Transcription Misc Initialization AutomaticDifferentation TestTools Models Transcription ${CMAKE_SOURCE_DIR})
target_link_libraries(iLQR_Include INTERFACE Eigen3::Eigen)
target_compile_options(iLQR_Include INTERFACE
  $<$<CXX_COMPILER_ID:MSVC>:/utf-8>
)
target_link_libraries(iLQR PUBLIC iLQR_Include)

if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
  # MSVC-like frontend (e.g. clang-cl).
  add_compile_options(/W4 /Gy)
  add_link_options(/OPT:REF /OPT:ICF)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  add_compile_options(-Wall -Wpedantic -Wextra -ffunction-sections -fdata-sections)
  add_link_options(-Wl,--gc-sections -Wl,--start-group -Wl,--end-group)
endif()

# 启用测试
enable_testing()

# 添加测试子目录
add_subdirectory(Tests)