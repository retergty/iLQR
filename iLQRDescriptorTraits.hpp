#pragma once

#include <array>

#include "Types.hpp"
#include "iLQRDescriptor.hpp"

/**
 * @brief Centralized view of an iLQR descriptor.
 *
 * Components should depend on this traits type when they need to read common
 * descriptor fields, instead of depending on the higher-level iLQRTypes
 * aggregation.
 */
template <typename Descriptor>
struct iLQRDescriptorTraits {
  using Descriptor_t = Descriptor;
  using Scalar = typename Descriptor::Scalar;
  using TranscriptionConfig = typename Descriptor::TranscriptionConfig;
  using Dims = typename Descriptor::Dims;
  using Horizon = typename Descriptor::Horizon;
  using ConstraintConfig = typename Descriptor::ConstraintConfig;

  using StateConstraintConfig = typename ConstraintConfig::State;
  using StateInputConstraintConfig = typename ConstraintConfig::StateInput;
  using FinalStateConstraintConfig = typename ConstraintConfig::FinalState;

  static constexpr int XDim = Dims::XDim;
  static constexpr int UDim = Dims::UDim;
  static constexpr std::size_t PredictLength = Horizon::PredictLength;

  static constexpr int StateEq = ConstraintConfig::StateEq;
  static constexpr int StateIneq = ConstraintConfig::StateIneq;
  static constexpr int StateInputEq = ConstraintConfig::StateInputEq;
  static constexpr int StateInputIneq = ConstraintConfig::StateInputIneq;
  static constexpr int FinalStateEq = ConstraintConfig::FinalStateEq;
  static constexpr int FinalStateIneq = ConstraintConfig::FinalStateIneq;

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

  using IntermediateStageConstraintLayout_t =
      IntermediateStageConstraintLayout<ConstraintConfig>;
  using FinalStageConstraintLayout_t =
      FinalStageConstraintLayout<ConstraintConfig>;
};
