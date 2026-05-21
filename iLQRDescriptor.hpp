#pragma once
#include <cstddef>

template <int X_, int U_>
struct Dimensions {
  static constexpr int XDim = X_;
  static constexpr int UDim = U_;
};

template <std::size_t N_>
struct Horizon {
  static constexpr std::size_t PredictLength = N_;
};

template <typename Dims_, typename Horizon_>
struct TranscriptionConfig {
  using Dims = Dims_;
  using Horizon = Horizon_;
  static constexpr int XDim = Dims::XDim;
  static constexpr int UDim = Dims::UDim;
  static constexpr std::size_t PredictLength = Horizon::PredictLength;
};

struct DefaultConstraintPolicy {};

template <int Eq_, int Ineq_>
struct ConstraintLayout {
  static constexpr int Eq = Eq_;
  static constexpr int Ineq = Ineq_;
};

template <typename Layout_, typename Policy_ = DefaultConstraintPolicy>
struct StateConstraintConfig {
  using Layout = Layout_;
  using Policy = Policy_;
};

template <typename Layout_, typename Policy_ = DefaultConstraintPolicy>
struct StateInputConstraintConfig {
  using Layout = Layout_;
  using Policy = Policy_;
};

template <typename Layout_, typename Policy_ = DefaultConstraintPolicy>
struct FinalStateConstraintConfig {
  using Layout = Layout_;
  using Policy = Policy_;
};

template <
    typename State_ = StateConstraintConfig<ConstraintLayout<0, 0>>,
    typename StateInput_ = StateInputConstraintConfig<ConstraintLayout<0, 0>>,
    typename FinalState_ = FinalStateConstraintConfig<ConstraintLayout<0, 0>>>
struct ConstraintConfig {
  using State = State_;
  using StateInput = StateInput_;
  using FinalState = FinalState_;

  static constexpr int StateEq = State::Layout::Eq;
  static constexpr int StateIneq = State::Layout::Ineq;
  static constexpr int StateInputEq = StateInput::Layout::Eq;
  static constexpr int StateInputIneq = StateInput::Layout::Ineq;
  static constexpr int FinalStateEq = FinalState::Layout::Eq;
  static constexpr int FinalStateIneq = FinalState::Layout::Ineq;
};

template <int StateEq_ = 0, int StateIneq_ = 0, int StateInputEq_ = 0,
          int StateInputIneq_ = 0>
struct StageConstraintLayout {
  static constexpr int StateEq = StateEq_;
  static constexpr int StateIneq = StateIneq_;
  static constexpr int StateInputEq = StateInputEq_;
  static constexpr int StateInputIneq = StateInputIneq_;
};

template <typename ConstraintConfig_>
using IntermediateStageConstraintLayout = StageConstraintLayout<
    ConstraintConfig_::StateEq, ConstraintConfig_::StateIneq,
    ConstraintConfig_::StateInputEq, ConstraintConfig_::StateInputIneq>;

template <typename ConstraintConfig_>
using FinalStageConstraintLayout =
    StageConstraintLayout<ConstraintConfig_::FinalStateEq,
                          ConstraintConfig_::FinalStateIneq, 0, 0>;

template <typename Scalar_, typename TranscriptionConfig_,
          typename ConstraintConfig_ = ConstraintConfig<>>
struct iLQRDescriptor {
  using Scalar = Scalar_;
  using TranscriptionConfig = TranscriptionConfig_;
  using Dims = typename TranscriptionConfig::Dims;
  using Horizon = typename TranscriptionConfig::Horizon;
  using ConstraintConfig = ConstraintConfig_;
};