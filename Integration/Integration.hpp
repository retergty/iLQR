#pragma once
#include "Types.hpp"
#include "Observer.hpp"
/**
 * The base class for autonomous system dynamics.
 */
template <typename Scalar, int XDimisions>
class OdeBase
{
public:
  /** Default constructor */
  OdeBase() = default;

  /** Default destructor */
  virtual ~OdeBase() = default;

  /**
   * Computes the autonomous system dynamics.
   * @param [in] t: Current time.
   * @param [in] x: Current state.
   * @return Current state time derivative
   */
  virtual Vector<Scalar, XDimisions> computeFlowMap(Scalar t, const Vector<Scalar, XDimisions> &x) const = 0;
};

/**
 * @brief The IntegratorType enum
 * Enum used in selecting a specific integrator.
 */
enum class IntegratorType
{
  EULER,
  ODE45,
  ODE45_OCS2,
  ADAMS_BASHFORTH,
  BULIRSCH_STOER,
  MODIFIED_MIDPOINT,
  RK4,
  RK5_VARIABLE,
  ADAMS_BASHFORTH_MOULTON
};

/**
 * The interface class for integration of autonomous systems.
 */
template <typename Scalar, int XDimisions>
class IntegratorBase
{
public:
  /**
   * Default constructor
   */
  IntegratorBase() = default;

  /**
   * Default destructor
   */
  virtual ~IntegratorBase() = default;

  /**
   * Equidistant integration based on initial and final time as well as step length.
   *
   * @param [in] system: System dynamics
   * @param [in] observer: Observer
   * @param [in] initialState: Initial state.
   * @param [in] startTime: Initial time.
   * @param [in] finalTime: Final time.
   * @param [in] dt: Time step.
   */
  virtual void integrateConst(OdeBase<Scalar, XDimisions> &system, Observer<Scalar, XDimisions> &observer, const Vector<Scalar, XDimisions> &initialState, const Scalar startTime, const Scalar finalTime, const Scalar dt) = 0;

  /**
   * Adaptive time integration based on start time and final time.
   *
   * @param [in] system: System dynamics
   * @param [in] observer: Observer
   * @param [in] initialState: Initial state.
   * @param [in] startTime: Initial time.
   * @param [in] finalTime: Final time.
   * @param [in] dtInitial: Initial time step.
   * @param [in] AbsTol: The absolute tolerance error for ode solver.
   * @param [in] RelTol: The relative tolerance error for ode solver.
   */
  virtual void integrateAdaptive(OdeBase<Scalar, XDimisions> &system, Observer<Scalar, XDimisions> &observer, const Vector<Scalar, XDimisions> &initialState, const Scalar startTime, const Scalar finalTime, const Scalar dtInitial = 0.01,
                                 const Scalar AbsTol = 1e-6, const Scalar RelTol = 1e-3) = 0;

  // /**
  //  * Output integration based on a given time trajectory.
  //  *
  //  * @param [in] system: System dynamics
  //  * @param [in] observer: Observer
  //  * @param [in] initialState: Initial state.
  //  * @param [in] beginTimeItr: The iterator to the beginning of the time stamp trajectory.
  //  * @param [in] endTimeItr: The iterator to the end of the time stamp trajectory.
  //  * @param [in] dtInitial: Initial time step.
  //  * @param [in] AbsTol: The absolute tolerance error for ode solver.
  //  * @param [in] RelTol: The relative tolerance error for ode solver.
  //  */
  // void integrateTimes(OdeBase &system, Observer &observer, const vector_t &initialState,
  //                     typename scalar_array_t::const_iterator beginTimeItr, typename scalar_array_t::const_iterator endTimeItr,
  //                     scalar_t dtInitial = 0.01, scalar_t AbsTol = 1e-6, scalar_t RelTol = 1e-3,
  //                     int maxNumSteps = std::numeric_limits<int>::max());
};