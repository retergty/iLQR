#pragma once
#include <array>
#include "LinearController.hpp"
/**
 * This class contains the primal problem's solution.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength>
struct PrimalSolution
{
  /** Constructor */
  PrimalSolution() = default;

  /** Destructor */
  ~PrimalSolution() = default;

  /** Copy constructor */
  PrimalSolution(const PrimalSolution &other)
      : timeTrajectory_(other.timeTrajectory_),
        stateTrajectory_(other.stateTrajectory_),
        inputTrajectory_(other.inputTrajectory_),
        controller_(other.controller_)
  {
  }

  /** Copy Assignment */
  PrimalSolution &operator=(const PrimalSolution &other)
  {
    timeTrajectory_ = other.timeTrajectory_;
    stateTrajectory_ = other.stateTrajectory_;
    inputTrajectory_ = other.inputTrajectory_;
    controller_ = other.controller_;
    return *this;
  }

  /** Move constructor */
  PrimalSolution(PrimalSolution &&other) noexcept = default;

  /** Move Assignment */
  PrimalSolution &operator=(PrimalSolution &&other) noexcept = default;

  /** Swap */
  void swap(PrimalSolution &other)
  {
    timeTrajectory_.swap(other.timeTrajectory_);
    stateTrajectory_.swap(other.stateTrajectory_);
    inputTrajectory_.swap(other.inputTrajectory_);
    controller_.swap(other.controller_);
  }

  void clear()
  {
    controller_.clear();
    for (size_t i = 0; i < PredictLength + 1; ++i)
    {
      timeTrajectory_[i] = 0;
      stateTrajectory_[i].setZero();
      inputTrajectory_[i].setZero();
    }
  }
  std::array<Scalar, PredictLength + 1> timeTrajectory_;
  std::array<Vector<Scalar, XDimisions>, PredictLength + 1> stateTrajectory_;
  std::array<Vector<Scalar, UDimisions>, PredictLength + 1> inputTrajectory_;
  LinearController<Scalar, XDimisions, UDimisions, PredictLength + 1> controller_;
};
