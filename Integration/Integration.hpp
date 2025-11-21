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
};