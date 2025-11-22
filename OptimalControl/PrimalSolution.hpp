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
        controllerPtr_(other.controllerPtr_)
  {
  }

  /** Copy Assignment */
  PrimalSolution &operator=(const PrimalSolution &other)
  {
    timeTrajectory_ = other.timeTrajectory_;
    stateTrajectory_ = other.stateTrajectory_;
    inputTrajectory_ = other.inputTrajectory_;
    controllerPtr_ = other.controllerPtr_;
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

    LinearController<Scalar, XDimisions, UDimisions, PredictLength + 1> *tmp = other.controllerPtr_;
    other.controllerPtr_ = controllerPtr_;
    controllerPtr_ = tmp;
  }

  std::array<Scalar, PredictLength + 1> timeTrajectory_;
  std::array<Vector<Scalar, XDimisions>, PredictLength + 1> stateTrajectory_;
  std::array<Vector<Scalar, UDimisions>, PredictLength + 1> inputTrajectory_;
  LinearController<Scalar, XDimisions, UDimisions, PredictLength + 1> *controllerPtr_{nullptr};
};
