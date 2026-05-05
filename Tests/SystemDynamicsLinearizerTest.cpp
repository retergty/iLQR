#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>

#include "LinearSystemDynamics.hpp"
#include "SystemDynamicsLinearizer.hpp"
#include "SystemDynamicsBase.hpp"
#include "Types.hpp"

namespace
{

constexpr double kTolerance = 1e-5;
constexpr double kEpsilon = 1e-10;
constexpr int kStateDim = 2;
constexpr int kInputDim = 1;
constexpr double kPi = 3.14159265358979323846;

using Scalar = double;
using StateVector = Vector<Scalar, kStateDim>;
using InputVector = Vector<Scalar, kInputDim>;
using StateMatrix = Matrix<Scalar, kStateDim, kStateDim>;
using InputMatrix = Matrix<Scalar, kStateDim, kInputDim>;
using LinearApproximation = VectorFunctionLinearApproximation<Scalar, kStateDim, kStateDim, kInputDim>;
using ControlledSystem = ControlledSystemBase<Scalar, kStateDim, kInputDim>;
using SystemDynamics = SystemDynamicsBase<Scalar, kStateDim, kInputDim>;
using LinearSystem = LinearSystemDynamics<Scalar, kStateDim, kInputDim>;

bool derivativeChecker(SystemDynamics &sys1, SystemDynamics &sys2, Scalar tolerance, Scalar t, const StateVector &x, const InputVector &u);

// Pendulum system, theta = 0 is upright.
class PendulumSystem final : public SystemDynamics
{
public:
    PendulumSystem() = default;
    ~PendulumSystem() override = default;

    StateVector computeFlowMap(Scalar t, const StateVector &x, const InputVector &u) const override
    {
        (void)t;
        StateVector dfdt;
        dfdt << x(1), std::sin(x(0)) + 0.1 * u(0);
        return dfdt;
    }

    LinearApproximation linearApproximation(Scalar t, const StateVector &x, const InputVector &u) override
    {
        (void)t;
        LinearApproximation linearDynamics;
        linearDynamics.f = computeFlowMap(t, x, u);
        linearDynamics.dfdx << 0.0, 1.0,
                               std::cos(x(0)), 0.0;
        linearDynamics.dfdu << 0.0,
                               0.1;
        return linearDynamics;
    }
};

TEST(testSystemDynamicsLinearizer, testDerivativeChecker)
{
    const Scalar time = 0.0;
    const StateVector state = StateVector::Zero();
    const InputVector input = InputVector::Zero();

    StateMatrix A;
    InputMatrix B;
    A << 0.6, 1.2,
         -0.8, 3.4;
    B << 1.0,
         1.0;

    LinearSystem linSys(A, B);

    A(0, 0) = 0.0;
    B(0, 0) = 0.0;
    LinearSystem alteredSys(A, B);

    ASSERT_FALSE(derivativeChecker(linSys, alteredSys, kTolerance, time, state, input));
}

TEST(testSystemDynamicsLinearizer, testLinearSystem)
{
    const Scalar time = 0.0;
    const StateVector state = StateVector::Zero();
    const InputVector input = InputVector::Zero();

    StateMatrix A;
    InputMatrix B;
    A << 0.6, 1.2,
         -0.8, 3.4;
    B << 1.0,
         1.0;

    LinearSystem linSys(A, B);

    SystemDynamicsLinearizer<Scalar, kStateDim, kInputDim> linearizedSys(
        std::unique_ptr<ControlledSystem>(std::make_unique<LinearSystem>(linSys)),
        /*doubleSidedDerivative=*/true,
        /*isSecondOrderSystem=*/false,
        kEpsilon);
    ASSERT_TRUE(derivativeChecker(linSys, linearizedSys, kTolerance, time, state, input));
}

TEST(testSystemDynamicsLinearizer, testPendulum)
{
    std::srand(0);

    constexpr std::size_t kDivisions = 1000;
    const Scalar maxDeg = 180.0;
    const Scalar toRads = kPi / 180.0;
    const Scalar t = 0.0;

    std::array<StateVector, kDivisions> testStates{};
    for (std::size_t i = 0; i < kDivisions; ++i)
    {
        const Scalar alpha = static_cast<Scalar>(i) / static_cast<Scalar>(kDivisions - 1);
        testStates[i] << alpha * toRads * maxDeg, 0.0;
    }

    const InputVector input = InputVector::Random();

    PendulumSystem nonLinSys;

    SystemDynamicsLinearizer<Scalar, kStateDim, kInputDim> linearizedSys(
        std::unique_ptr<ControlledSystem>(std::make_unique<PendulumSystem>(nonLinSys)),
        /*doubleSidedDerivative=*/true,
        /*isSecondOrderSystem=*/false,
        kEpsilon);
    for (const auto &state : testStates)
    {
        ASSERT_TRUE(derivativeChecker(nonLinSys, linearizedSys, kTolerance, t, state, input));
    }
}

bool derivativeChecker(SystemDynamics &sys1, SystemDynamics &sys2, Scalar tolerance, Scalar t, const StateVector &x, const InputVector &u)
{
    const auto derivatives1 = sys1.linearApproximation(t, x, u);
    const auto derivatives2 = sys2.linearApproximation(t, x, u);
    const Scalar AError = (derivatives1.dfdx - derivatives2.dfdx).lpNorm<Eigen::Infinity>();
    const Scalar BError = (derivatives1.dfdu - derivatives2.dfdu).lpNorm<Eigen::Infinity>();
    return tolerance > std::max(AError, BError);
}

} // namespace