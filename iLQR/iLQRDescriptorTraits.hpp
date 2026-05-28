#pragma once

#include <array>

#include "OptimalControlData/Reference.hpp"
#include "Types.hpp"
#include "iLQR/iLQRDescriptor.hpp"

/**
 * @brief iLQR descriptor 的集中视图。
 *
 * 组件需要读取公共 descriptor 字段时，应依赖该 traits 类型，
 * 而不是依赖更高层的 iLQRTypes
 * 聚合。
 */
template <typename Descriptor>
struct iLQRDescriptorTraits {
  using Descriptor_t = Descriptor;
  using Scalar = typename Descriptor::Scalar;
  using TranscriptionConfig = typename Descriptor::TranscriptionConfig;
  using Dims = typename Descriptor::Dims;
  using Horizon = typename Descriptor::Horizon;
  using DynamicsMode = typename TranscriptionConfig::DynamicsMode;
  using ConstraintConfig = typename Descriptor::ConstraintConfig;

  using StateConstraintConfig = typename ConstraintConfig::State;
  using StateInputConstraintConfig = typename ConstraintConfig::StateInput;
  using FinalStateConstraintConfig = typename ConstraintConfig::FinalState;
  using StateEqLayout = typename ConstraintConfig::StateEqLayout;
  using StateIneqLayout = typename ConstraintConfig::StateIneqLayout;
  using StateInputEqLayout = typename ConstraintConfig::StateInputEqLayout;
  using StateInputIneqLayout = typename ConstraintConfig::StateInputIneqLayout;
  using FinalStateEqLayout = typename ConstraintConfig::FinalStateEqLayout;
  using FinalStateIneqLayout = typename ConstraintConfig::FinalStateIneqLayout;

  static constexpr int XDim = Dims::XDim;
  static constexpr int UDim = Dims::UDim;
  static constexpr std::size_t PredictLength = Horizon::PredictLength;

  static constexpr int StateEqDim = StateEqLayout::TotalDim;
  static constexpr int StateIneqDim = StateIneqLayout::TotalDim;
  static constexpr int StateInputEqDim = StateInputEqLayout::TotalDim;
  static constexpr int StateInputIneqDim = StateInputIneqLayout::TotalDim;
  static constexpr int FinalStateEqDim = FinalStateEqLayout::TotalDim;
  static constexpr int FinalStateIneqDim = FinalStateIneqLayout::TotalDim;

  using StateVector_t = Vector<Scalar, XDim>;
  using InputVector_t = Vector<Scalar, UDim>;
  using StateMatrix_t = Matrix<Scalar, XDim, XDim>;
  using InputMatrix_t = Matrix<Scalar, UDim, UDim>;
  using InputStateMatrix_t = Matrix<Scalar, UDim, XDim>;
  using StateInputMatrix_t = Matrix<Scalar, XDim, UDim>;
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
  using TargetTrajectories_t = TargetTrajectories<Scalar, TranscriptionConfig>;

  using IntermediateStageConstraintLayout_t =
      IntermediateStageConstraintLayout<ConstraintConfig>;
  using FinalStageConstraintLayout_t =
      FinalStageConstraintLayout<ConstraintConfig>;
};
