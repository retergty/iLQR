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
  /** This value determines the absolute tolerance error for ode solvers. */
  Scalar absTolODE = 1e-9;
  /** This value determines the relative tolerance error for ode solvers. */
  Scalar relTolODE = 1e-6;
  /** This value determines the maximum number of integration points per a second for ode solvers. */
  Scalar maxNumStepsPerSecond = 10000;
  /** The integration time step used in the fixed time-step rollout methods */
  Scalar timeStep = 1e-2;
  /** Rollout integration scheme type */
  IntegratorType integratorType = IntegratorType::ODE45;

  /** Whether to run controller again after integration to construct input trajectory */
  bool reconstructInputTrajectory = true;

  /** Which of the RootFinding algorithms to use in StateRollout
   * 		0:		Anderson & Björck		(default)
   * 		1:		Pegasus
   * 		2:		Illinois
   * 		3:		Regula Falsi
   */
  RootFinderType rootFindingAlgorithm = RootFinderType::ANDERSON_BJORCK;
  /** Whether to use the trajectory spreading controller in state triggered rollout */
  bool useTrajectorySpreadingController = false;
};

/**
 * This class is an interface class for forward rollout of the system dynamics.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t ArrayLen>
class RolloutBase
{
public:
  /**
   * Default constructor.
   *
   * @param [in] rolloutSettings: The rollout settings.
   */
  explicit RolloutBase(const RolloutSettings<Scalar> &rolloutSettings) : rolloutSettings_(rolloutSettings) {}

  /**
   * Default destructor.
   */
  virtual ~RolloutBase() = default;

  /**
   * Returns the rollout settings.
   *
   * @return The rollout settings.
   */
  const RolloutSettings<Scalar> &settings() const { return rolloutSettings_; }

  /**
   * The kills the integrator inside the rollout.
   */
  virtual void abortRollout() {}

  /**
   * The enables the integrator inside the rollout to start again.
   */
  virtual void reactivateRollout() {}

  /**
   * Resets the simulator state to the initial state in the next runImpl.
   * @note This is relevant if a physics engine (e.g. RaiSim) is used.
   */
  virtual void resetRollout() {}

  /**
   * Forward integrate the system dynamics with given controller. It uses the given control policies and initial state,
   * to integrate the system dynamics in time period [initTime, finalTime].
   *
   * @param [in] initTime: The initial time.
   * @param [in] initState: The initial state.
   * @param [in] finalTime: The final time.
   * @param [in] controller: control policy.
   * @param [out] timeTrajectory: The time trajectory stamp.
   * @param [out] stateTrajectory: The state trajectory.
   * @param [out] inputTrajectory: The control input trajectory.
   *
   * @return numbers of array member
   */
  virtual int run(const Scalar initTime, const Vector<Scalar, XDimisions> &initState, const Scalar finalTime, ControllerBase<Scalar,XDimisions,UDimisions> *controller,
                                         std::array<Scalar, ArrayLen> &timeTrajectory, std::array<Vector<Scalar, XDimisions>, ArrayLen> &stateTrajectory, std::array<Vector<Scalar, UDimisions>, ArrayLen> &inputTrajectory) = 0;

protected:
  RolloutSettings<Scalar> rolloutSettings_;
};