#pragma once
#include "Types.hpp"
#include "Integration.hpp"
#include "LinearInterpolation.hpp"
#include <array>
#include <algorithm>
/**
 * This class is an interface class for the user defined target trajectories.
 */
template <typename Scalar, int XDimisions, int UDimisions, int ArrayLength>
struct TargetTrajectories
{
  TargetTrajectories() = default;

  TargetTrajectories(
    const std::array<Scalar, ArrayLength>& desiredTimeTrajectory,
    const std::array<Vector<Scalar, XDimisions>, ArrayLength>& desiredStateTrajectory,
    const std::array<Vector<Scalar, UDimisions>, ArrayLength>& desiredInputTrajectory)
    : timeTrajectory(desiredTimeTrajectory), stateTrajectory(desiredStateTrajectory), inputTrajectory(desiredInputTrajectory) {
  }

  TargetTrajectories(
    const std::array<Scalar, ArrayLength>& desiredTimeTrajectory,
    const std::array<Vector<Scalar, XDimisions>, ArrayLength>& desiredStateTrajectory)
    : timeTrajectory(desiredTimeTrajectory), stateTrajectory(desiredStateTrajectory)
  {
    for (auto& vec : inputTrajectory)
    {
      vec.setZero();
    }
  }

  void setTrajectory(const std::array<Scalar, ArrayLength>& desiredTimeTrajectory,
    const std::array<Vector<Scalar, XDimisions>, ArrayLength>& desiredStateTrajectory,
    const std::array<Vector<Scalar, UDimisions>, ArrayLength>& desiredInputTrajectory)
  {
    timeTrajectory = desiredTimeTrajectory;
    stateTrajectory = desiredStateTrajectory;
    inputTrajectory = desiredInputTrajectory;
  }

  void setTrajectory(const std::array<Scalar, ArrayLength>& desiredTimeTrajectory,
    const std::array<Vector<Scalar, XDimisions>, ArrayLength>& desiredStateTrajectory)
  {
    timeTrajectory = desiredTimeTrajectory;
    stateTrajectory = desiredStateTrajectory;
    for (auto& vec : inputTrajectory)
    {
      vec.setZero();
    }
  }

  Vector<Scalar, XDimisions> getDesiredState(const int index) const
  {
    assert(index >= 0 && index < ArrayLength);
    return stateTrajectory[index];
  }
  template <int Index>
  Vector<Scalar, XDimisions> getDesiredState() const
  {
    assert(Index >= 0 && Index < ArrayLength);
    return stateTrajectory[Index];
  }
  Vector<Scalar, XDimisions> getDesiredState(const Scalar time) const
  {
    return LinearInterpolation::interpolate(time, timeTrajectory, stateTrajectory);
  }

  Vector<Scalar, UDimisions> getDesiredInput(const int index) const
  {
    assert(index >= 0 && index < ArrayLength);
    return inputTrajectory[index];
  }

  template <int Index>
  Vector<Scalar, UDimisions> getDesiredInput() const
  {
    assert(Index >= 0 && Index < ArrayLength);
    return inputTrajectory[Index];
  }
  Vector<Scalar, UDimisions> getDesiredInput(const Scalar time) const
  {
    return LinearInterpolation::interpolate(time, timeTrajectory, inputTrajectory);
  }

  std::array<Scalar, ArrayLength> timeTrajectory;
  std::array<Vector<Scalar, XDimisions>, ArrayLength> stateTrajectory;
  std::array<Vector<Scalar, UDimisions>, ArrayLength> inputTrajectory;
};