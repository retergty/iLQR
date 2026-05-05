#pragma once

#include <memory>

#include "ControlledSystemBase.hpp"
#include "SystemDynamicsBase.hpp"
#include "FiniteDifferenceMethods.hpp"
/**
 * A class for linearizing system dynamics. The linearized system dynamics is defined as: \n
 *
 * - Linearized system:   \f$ dx/dt = A(t) \delta x + B(t) \delta u \f$ \n
 */
template <typename Scalar, int XDimisions, int UDimisions>
class SystemDynamicsLinearizer final : public SystemDynamicsBase<Scalar, XDimisions, UDimisions> {
public:
  using StateVector_t = Vector<Scalar, XDimisions>;
  using InputVector_t = Vector<Scalar, UDimisions>;
  using VectorFunctionLinearApproximation_t = VectorFunctionLinearApproximation<Scalar, XDimisions, XDimisions, UDimisions>;
  /** Constructor */
  explicit SystemDynamicsLinearizer(ControlledSystemBase* nonlinearSystemPtr, bool doubleSidedDerivative = true,
    bool isSecondOrderSystem = false, Scalar eps = Eigen::NumTraits<Scalar>::epsilon()) : SystemDynamicsBase<Scalar, XDimisions, UDimisions>(),
    controlledSystemPtr_(nonlinearSystemPtr),
    doubleSidedDerivative_(doubleSidedDerivative),
    isSecondOrderSystem_(isSecondOrderSystem),
    eps_(eps)
  {
  }

  ~SystemDynamicsLinearizer() override = default;

  StateVector_t computeFlowMap(Scalar time, const StateVector_t& state, const InputVector_t& input) override const
  {
    return controlledSystemPtr_->computeFlowMap(time, state, input);
  }

  VectorFunctionLinearApproximation_t linearApproximation(Scalar t, const StateVector_t& x, const InputVector_t& u) override
  {
    VectorFunctionLinearApproximation_t linearDynamics;
    linearDynamics.f = controlledSystemPtr_->computeFlowMap(t, x, u);
    linearDynamics.dfdx = finiteDifferenceDerivativeState(*controlledSystemPtr_, t, x, u, eps_, doubleSidedDerivative_, isSecondOrderSystem_);
    linearDynamics.dfdu = finiteDifferenceDerivativeInput(*controlledSystemPtr_, t, x, u, eps_, doubleSidedDerivative_, isSecondOrderSystem_);
    return linearDynamics;
  }

private:
  SystemDynamicsLinearizer(const SystemDynamicsLinearizer& other);

  ControlledSystemBase<Scalar, XDimisions, UDimisions>* controlledSystemPtr_;
  bool doubleSidedDerivative_;
  bool isSecondOrderSystem_;
  Scalar eps_;
};

