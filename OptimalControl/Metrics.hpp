#pragma once

#include "Types.hpp"
#include "LagrangianMetrics.hpp"
#include "Numerics.hpp"
#include <cmath>
/**
 * The collection of cost, dynamics violation, constraints, and LagrangianMetrics structure for all possible constraint
 * terms (handled by Lagrangian method) in a point of time.
 *     cost : The total cost in a particular time point.
 *     dynamicsViolation : The vector of dynamics violation.
 *     stateEqLagrangian : An array of state equality constraint terms handled by Lagrangian method.
 *     stateIneqLagrangian : An array of state inequality constraint terms handled by Lagrangian method.
 *     stateInputEqLagrangian : An array of state-input equality constraint terms handled by Lagrangian method.
 *     stateInputIneqLagrangian : An array of state-input inequality constraint terms handled by Lagrangian method.
 */
template <typename Scalar, int XDimisions, int UDimisions, int StateEqConstrains, int StateIneqConstrains, int StateInputEqConstrains, int StateInputIneqConstrains>
struct Metrics
{
    // Cost
    Scalar cost;

    // Dynamics violation
    Vector<Scalar, XDimisions> dynamicsViolation;

    // Lagrangians
    std::array<LagrangianMetrics<Scalar>, StateEqConstrains> stateEqLagrangian;
    std::array<LagrangianMetrics<Scalar>, StateIneqConstrains> stateIneqLagrangian;
    std::array<LagrangianMetrics<Scalar>, StateInputEqConstrains> stateInputEqLagrangian;
    std::array<LagrangianMetrics<Scalar>, StateInputIneqConstrains> stateInputIneqLagrangian;

    /** Exchanges the values of Metrics */
    void swap(Metrics& other)
    {
        // Cost
        std::swap(cost, other.cost);
        // Dynamics violation
        dynamicsViolation.swap(other.dynamicsViolation);
        // Lagrangians
        stateEqLagrangian.swap(other.stateEqLagrangian);
        stateIneqLagrangian.swap(other.stateIneqLagrangian);
        stateInputEqLagrangian.swap(other.stateInputEqLagrangian);
        stateInputIneqLagrangian.swap(other.stateInputIneqLagrangian);
    }

    // /** Clears the value of the Metrics */
    // void clear();

    // /** Returns true if *this is approximately equal to other, within the precision determined by prec. */
    // bool isApprox(const Metrics &other, Scalar prec = 1e-8) const
    // {

    // }
};

/** Sums penalties of an array of LagrangianMetrics. */
template<typename Scalar, size_t ArrayLen>
inline Scalar sumPenalties(const std::array<LagrangianMetrics<Scalar>, ArrayLen>& metricsArray)
{
    Scalar s = 0.0;
    //std::for_each(metricsArray.begin(), metricsArray.end(), [&s](const LagrangianMetrics<Scalar>& m) { s += m.penalty; });
    for (size_t i = 0;i < ArrayLen;++i)
    {
        s += metricsArray[i].penalty;
    }
    return s;
}

/** Computes the sum of squared norm of constraints of an array of LagrangianMetrics. */
template<typename Scalar, size_t ArrayLen>
inline Scalar constraintsSquaredNorm(const std::array<LagrangianMetrics<Scalar>, ArrayLen>& metricsArray) {
    Scalar s = 0.0;
    //std::for_each(metricsArray.begin(), metricsArray.end(), [&s](const LagrangianMetrics& m) { s += m.constraint.squaredNorm(); });
    for (int i = 0;i < ArrayLen;++i)
    {
        s += std::abs(metricsArray[i].constraint);
    }
    return s;
}

/** Computes the sum of squared norm of a vector of equality constraints violation. */
template<typename Scalar, int Dimisions>
inline Scalar getEqConstraintsSSE(const Vector<Scalar, Dimisions>& eqConstraint)
{
    static_assert(Dimisions >= 0);

    if constexpr (Dimisions == 0)
    {
        return 0;
    }
    else {
        return eqConstraint.squaredNorm();
    }
}