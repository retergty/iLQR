#pragma once
#include <cstddef>

template <int X_, int U_, std::size_t N_>
struct Dimensions {
  static constexpr int XDim = X_;
  static constexpr int UDim = U_;
  static constexpr std::size_t PredictLength = N_;
};

template <int Eq_ = 0, int Ineq_ = 0>
struct StateConstraints {
  static constexpr int EqCount = Eq_;
  static constexpr int IneqCount = Ineq_;
};

template <int Eq_ = 0, int Ineq_ = 0>
struct StateInputConstraints {
  static constexpr int EqCount = Eq_;
  static constexpr int IneqCount = Ineq_;
};

template <int Eq_ = 0, int Ineq_ = 0>
struct FinalStateConstraints {
  static constexpr int EqCount = Eq_;
  static constexpr int IneqCount = Ineq_;
};

template <typename State_ = StateConstraints<>,
          typename StateInput_ = StateInputConstraints<>,
          typename FinalState_ = FinalStateConstraints<>>
struct Constraints {
  using State_t = State_;
  using StateInput_t = StateInput_;
  using FinalState_t = FinalState_;

  static constexpr int StateEq = State_::EqCount;
  static constexpr int StateIneq = State_::IneqCount;
  static constexpr int StateInputEq = StateInput_::EqCount;
  static constexpr int StateInputIneq = StateInput_::IneqCount;
  static constexpr int FinalStateEq = FinalState_::EqCount;
  static constexpr int FinalStateIneq = FinalState_::IneqCount;
};

template <typename Scalar_, typename Dims_,
          typename ConstraintConfig_ = Constraints<>>
struct iLQRDescriptor {
  using Scalar = Scalar_;
  using Dims = Dims_;
  using Constraints = ConstraintConfig_;
};