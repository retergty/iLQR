#pragma once

#include "Types.hpp"
#include "Controller.hpp"
#include "Integration.hpp"

enum class RootFinderType
{
  ANDERSON_BJORCK,
  PEGASUS,
  ILLINOIS,
  REGULA_FALSI
};

/**
 * This structure contains the settings for forward rollout algorithms.
 */
template <typename Scalar>
struct RolloutSettings
{
  /** The integration time step used in the fixed time-step rollout methods */
  Scalar timeStep = 1e-2;

  /** Whether to run controller again after integration to construct input trajectory */
  bool reconstructInputTrajectory = true;
};

template <typename Scalar, int XDimisions, int UDimisions>
struct RolloutTrajectoryPointer
{
  RolloutTrajectoryPointer(Scalar* time_trajectory, Vector<Scalar, XDimisions>* state_trajectory, Vector<Scalar, UDimisions>* input_trajectory, int max_length)
    : timeTrajectory(time_trajectory), stateTrajectory(state_trajectory), inputTrajectory(input_trajectory), maxLength(max_length)
  {
  };
  Scalar* timeTrajectory;
  Vector<Scalar, XDimisions>* stateTrajectory;
  Vector<Scalar, UDimisions>* inputTrajectory;
  size_t maxLength;
};

/**
 * This class is an interface class for forward rollout of the system dynamics.
 */
template <typename Scalar, int XDimisions, int UDimisions>
class RolloutBase
{
public:
  using RolloutTrajectoryPointer_t = RolloutTrajectoryPointer<Scalar, XDimisions, UDimisions>;
  /**
   * Default constructor.
   *
   * @param [in] rolloutSettings: The rollout settings.
   */
  explicit RolloutBase() {}

  /**
   * Default destructor.
   */
  virtual ~RolloutBase() = default;

  /**
   * Returns the rollout settings.
   *
   * @return The rollout settings.
   */
  RolloutSettings<Scalar>& settings() { return rolloutSettings_; }

  /**
   * Forward integrate the system dynamics with given controller. It uses the given control policies and initial state,
   * to integrate the system dynamics in time period [initTime, finalTime].
   *
   * @param [in] initTime: The initial time.
   * @param [in] initState: The initial state.
   * @param [in] finalTime: The final time.
   * @param [in] controller: control policy.
   * @param [in out] RolloutTrajectoryPointer: The Rollout Trajectory.
   * @return numbers of array member
   */
  virtual int run(const Scalar initTime, const Vector<Scalar, XDimisions>& initState, const Scalar finalTime, ControllerBase<Scalar, XDimisions, UDimisions>* controller,
    RolloutTrajectoryPointer_t& trajectory) = 0;

  // /**
  //  * Forward integrate the system dynamics with given controller. It uses the given control policies and initial state,
  //  * to integrate the system dynamics in time period [initTime, finalTime].
  //  *
  //  * @param [in] initTime: The initial time.
  //  * @param [in] initState: The initial state.
  //  * @param [in] steps: numbers of steps.
  //  * @param [in] controller: control policy.
  //  * @param [out] timeTrajectory: The time trajectory stamp.
  //  * @param [out] stateTrajectory: The state trajectory.
  //  * @param [out] inputTrajectory: The control input trajectory.
  //  * @return the final time
  //  */
  // virtual int run(const Scalar initTime, const Vector<Scalar, XDimisions>& initState, const int steps, ControllerBase<Scalar, XDimisions, UDimisions>& controller,
  //   std::array<Scalar, ArrayLen>& timeTrajectory, std::array<Vector<Scalar, XDimisions>, ArrayLen>& stateTrajectory, std::array<Vector<Scalar, UDimisions>, ArrayLen>& inputTrajectory) = 0;

protected:
  RolloutSettings<Scalar> rolloutSettings_{};
};