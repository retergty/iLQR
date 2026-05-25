#pragma once
#include <cstddef>
#include <tuple>

template <int X_, int U_>
struct Dimensions {
  static constexpr int XDim = X_;
  static constexpr int UDim = U_;
};

template <std::size_t N_>
struct Horizon {
  static constexpr std::size_t PredictLength = N_;
};

/** @brief 连续时间动力学模式：模型提供 \f$\dot{x}=f(t,x,u)\f$。 */
struct ContinuousDynamics {};

/** @brief 离散时间动力学模式：模型提供 \f$x_{k+1}=f_d(t,x,u,dt)\f$。 */
struct DiscreteDynamics {};

template <typename Dims_, typename Horizon_,
          typename DynamicsMode_ = ContinuousDynamics>
struct TranscriptionConfig {
  using Dims = Dims_;
  using Horizon = Horizon_;
  using DynamicsMode = DynamicsMode_;
  static constexpr int XDim = Dims::XDim;
  static constexpr int UDim = Dims::UDim;
  static constexpr std::size_t PredictLength = Horizon::PredictLength;
};

struct DefaultConstraintPolicy {};

template <int CDim_>
struct ConstraintTerm {
  static constexpr int CDim = CDim_;
};

template <typename... Terms_>
struct ConstraintGroupLayout {
  using Terms = std::tuple<Terms_...>;
  static constexpr int NumTerms = sizeof...(Terms_);
  static constexpr int TotalDim = (Terms_::CDim + ... + 0);
  template <std::size_t I>
  using Term = std::tuple_element_t<I, Terms>;
};

template <typename Eq_ = ConstraintGroupLayout<>,
          typename Ineq_ = ConstraintGroupLayout<>>
struct ConstraintLayout {
  using Eq = Eq_;
  using Ineq = Ineq_;
  static constexpr int EqNumTerms = Eq::NumTerms;
  static constexpr int IneqNumTerms = Ineq::NumTerms;
  static constexpr int EqTotalDim = Eq::TotalDim;
  static constexpr int IneqTotalDim = Ineq::TotalDim;
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

template <typename State_ = StateConstraintConfig<ConstraintLayout<>>,
          typename StateInput_ = StateInputConstraintConfig<ConstraintLayout<>>,
          typename FinalState_ = FinalStateConstraintConfig<ConstraintLayout<>>>
struct ConstraintConfig {
  using State = State_;
  using StateInput = StateInput_;
  using FinalState = FinalState_;

  using StateEqLayout = typename State::Layout::Eq;
  using StateIneqLayout = typename State::Layout::Ineq;
  using StateInputEqLayout = typename StateInput::Layout::Eq;
  using StateInputIneqLayout = typename StateInput::Layout::Ineq;
  using FinalStateEqLayout = typename FinalState::Layout::Eq;
  using FinalStateIneqLayout = typename FinalState::Layout::Ineq;

  static constexpr int StateEq = State::Layout::EqTotalDim;
  static constexpr int StateIneq = State::Layout::IneqTotalDim;
  static constexpr int StateInputEq = StateInput::Layout::EqTotalDim;
  static constexpr int StateInputIneq = StateInput::Layout::IneqTotalDim;
  static constexpr int FinalStateEq = FinalState::Layout::EqTotalDim;
  static constexpr int FinalStateIneq = FinalState::Layout::IneqTotalDim;
};

template <typename StateEq_ = ConstraintGroupLayout<>,
          typename StateIneq_ = ConstraintGroupLayout<>,
          typename StateInputEq_ = ConstraintGroupLayout<>,
          typename StateInputIneq_ = ConstraintGroupLayout<>>
struct StageConstraintLayout {
  using StateEqLayout = StateEq_;
  using StateIneqLayout = StateIneq_;
  using StateInputEqLayout = StateInputEq_;
  using StateInputIneqLayout = StateInputIneq_;

  static constexpr int StateEq = StateEqLayout::TotalDim;
  static constexpr int StateIneq = StateIneqLayout::TotalDim;
  static constexpr int StateInputEq = StateInputEqLayout::TotalDim;
  static constexpr int StateInputIneq = StateInputIneqLayout::TotalDim;

  static constexpr int StateEqNumTerms = StateEqLayout::NumTerms;
  static constexpr int StateIneqNumTerms = StateIneqLayout::NumTerms;
  static constexpr int StateInputEqNumTerms = StateInputEqLayout::NumTerms;
  static constexpr int StateInputIneqNumTerms = StateInputIneqLayout::NumTerms;
};

template <typename ConstraintConfig_>
using IntermediateStageConstraintLayout =
    StageConstraintLayout<typename ConstraintConfig_::StateEqLayout,
                          typename ConstraintConfig_::StateIneqLayout,
                          typename ConstraintConfig_::StateInputEqLayout,
                          typename ConstraintConfig_::StateInputIneqLayout>;

template <typename ConstraintConfig_>
using FinalStageConstraintLayout =
    StageConstraintLayout<typename ConstraintConfig_::FinalStateEqLayout,
                          typename ConstraintConfig_::FinalStateIneqLayout,
                          ConstraintGroupLayout<>, ConstraintGroupLayout<>>;

template <typename Scalar_, typename TranscriptionConfig_,
          typename ConstraintConfig_ = ConstraintConfig<>>
struct iLQRDescriptor {
  using Scalar = Scalar_;
  using TranscriptionConfig = TranscriptionConfig_;
  using Dims = typename TranscriptionConfig::Dims;
  using Horizon = typename TranscriptionConfig::Horizon;
  using ConstraintConfig = ConstraintConfig_;
};