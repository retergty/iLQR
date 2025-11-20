#pragma once
#include "Dynamics.hpp"
/**
 * The system dynamics and linearization class.
 * The linearized system flow map is defined as: \n
 * \f$ dx/dt = A(t) \delta x + B(t) \delta u \f$ \n
 * The linearized system jump map is defined as: \n
 * \f$ x^+ = G \delta x + H \delta u \f$ \n
 */
template <typename Scalar, int XDimisions, int UDimisions>
class SystemDynamicsBase : public ControlledSystemBase<Scalar, XDimisions, UDimisions>
{
public:
  /**
   * Constructor
   *
   */
  SystemDynamicsBase() = default;

  /** Default destructor */
  ~SystemDynamicsBase() override = default;

  /**
   * Computes the linear approximation.
   *
   * @param [in] t: The current time.
   * @param [in] x: The current state.
   * @param [in] u: The current input.
   * @param [in] preComp: pre-computation module, safely ignore this parameter if not used.
   *                      @see PreComputation class documentation.
   * @return The state time derivative linear approximation.
   */
  virtual VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>
  linearApproximation(Scalar t, const Vector<Scalar, XDimisions> &x, const Vector<Scalar, UDimisions> &u) = 0;
};
