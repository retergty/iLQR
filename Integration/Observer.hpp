#pragma once

#include "Types.hpp"

/**
 * The Observer class stores data in given containers.
 */
template <typename Scalar, int XDimisions>
class Observer
{
public:
  /**
   * Constructor.
   *
   * @param stateTrajectoryPtr: A pinter to an state trajectory container to store resulting state trajectory.
   * @param timeTrajectoryPtr: A pinter to an time trajectory container to store resulting time trajectory.
   */
  explicit Observer(int length, Vector<Scalar, XDimisions>* stateTrajectoryPtr = nullptr, Scalar* timeTrajectoryPtr = nullptr) : stateTrajectoryPtr_(stateTrajectoryPtr), timeTrajectoryPtr_(timeTrajectoryPtr), length_(length) {};

  /**
   * Default destructor.
   */
  ~Observer() = default;

  /**
   * Observe function to retrieve the variable of interest.
   * @param [in] state: Current state.
   * @param [in] time: Current time.
   */
  void observe(const std::array<Scalar, XDimisions>& state, const Scalar time)
  {
    if (now_ >= length_)
      return;
    if (timeTrajectoryPtr_ != nullptr)
      timeTrajectoryPtr_[now_] = time;
    if (stateTrajectoryPtr_ != nullptr)
      stateTrajectoryPtr_[now_] = state;
    now_++;
    return;
  }

  void clear()
  {
    now_ = 0;
  }

  int getCount() const
  {
    return now_;
  }
private:
  Vector<Scalar, XDimisions>* timeTrajectoryPtr_;
  Scalar* stateTrajectoryPtr_;
  int length_;
  int now_{ 0 };
};
