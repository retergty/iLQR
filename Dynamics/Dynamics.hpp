#pragma once
#include "Types.hpp"
#include "Integration.hpp"
#include "Controller.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"

/**
 * The base class for non-autonomous system dynamics.
 */
template <typename Scalar, int XDimisions, int UDimisions>
class ControlledSystemBase : public OdeBase<Scalar, XDimisions>
{
public:
  /**
   * Constructor
   *
   */
  ControlledSystemBase() = default;

  /** Default destructor */
  ~ControlledSystemBase() override = default;

  /**
   * Computes the flow map of a system.
   *
   * @param [in] t: The current time.
   * @param [in] x: The current state.
   * @return The state time derivative.
   */
  Vector<Scalar, XDimisions> computeFlowMap(Scalar t, const Vector<Scalar, XDimisions> &x) const override final
  {
    assert(controllerPtr_ != nullptr);
    const Vector<Scalar, UDimisions> u = controllerPtr_->computeInput(t, x);
    return computeFlowMap(t, x, u);
  }

  /**
   * Computes the flow map of a system with exogenous input.
   *
   * @param [in] t: The current time.
   * @param [in] x: The current state.
   * @param [in] u: The current input.
   * @return The state time derivative.
   */
  virtual Vector<Scalar, XDimisions> computeFlowMap(Scalar t, const Vector<Scalar, XDimisions> &x, const Vector<Scalar, UDimisions> &u) const = 0;

  void setController(ControllerBase<Scalar, XDimisions, UDimisions> *controllerPtr)
  {
    controllerPtr_ = controllerPtr;
  }
  const ControllerBase<Scalar, XDimisions, UDimisions> *controllerPtr() const
  {
    return controllerPtr_;
  }

private:
  ControllerBase<Scalar, XDimisions, UDimisions> *controllerPtr_ = nullptr; //! pointer to controller
};
