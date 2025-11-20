#pragma once
#include "Cost.hpp"
#include <memory>
#include "IntrusiveList.hpp"

/**
 * Cost function combining a collection of cost terms.
 *
 * This class collects a variable number of cost terms and provides methods to get the
 * summed cost values and quadratic approximations. Each cost term can be accessed through its
 * string name and can be activated or deactivated.
 */
template <typename Scalar, int XDimisions, int ArrayLength>
class StateCostCollection
{
public:
  StateCostCollection() = default;
  virtual ~StateCostCollection() = default;

  /** Get state-only cost value */
  Scalar getValue(Scalar time, const Vector<Scalar, XDimisions> &state, const std::array<Scalar, ArrayLength> &timeTrajectories, const std::array<Vector<Scalar, XDimisions>, ArrayLength> &stateTrajectoies) const
  {
    Scalar cost = 0;
    for (auto it = list_.begin(); it != list_.end(); ++it)
    {
      cost += it->getValue(time, state, timeTrajectories, stateTrajectoies);
    }
    return cost;
  }

  /** Get state-only cost quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
  getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions> &state, const std::array<Scalar, ArrayLength> &timeTrajectories, const std::array<Vector<Scalar, XDimisions>, ArrayLength> &stateTrajectoies) const
  {
    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> cost_appro;
    for (auto it = list_.begin(); it != list_.end(); ++it)
    {
      cost_appro += it->getQuadraticApproximation(time, state, timeTrajectories, stateTrajectoies);
    }
    return cost_appro;
  }

  // add cost to list end
  void add(const StateCost<Scalar, XDimisions, ArrayLength> &cost)
  {
    list_.insert(list_.end(), cost);
  }

private:
  IntrusiveList<StateCost<Scalar, XDimisions, ArrayLength>> list_;
};

/**
 * State Input Cost function combining a collection of cost terms.
 *
 * This class collects a variable number of cost terms and provides methods to get the
 * summed cost values and quadratic approximations. Each cost term can be accessed through its
 * string name and can be activated or deactivated.
 */
template <typename Scalar, int XDimisions, int UDimisions, int ArrayLength>
class StateInputCostCollection
{
public:
  StateInputCostCollection() = default;
  ~StateInputCostCollection() = default;

  /** Get state-input cost value */
  Scalar getValue(
      Scalar time, const Vector<Scalar, XDimisions> &state, const Vector<Scalar, UDimisions> &input,
      const std::array<Scalar, ArrayLength> &timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength> &stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength> &inputTrajectory) const
  {
    Scalar cost = 0;
    for (auto it = list_.begin(); it != list_.end(); ++it)
    {
      cost += it->getValue(time, state, input, timeTrajectory, stateTrajectoy, inputTrajectory);
    }
    return cost;
  }

  /** Get state-input cost quadratic approximation */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>
  getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions> &state, const Vector<Scalar, UDimisions> &input,
                            const std::array<Scalar, ArrayLength> &timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength> &stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength> &inputTrajectory) const
  {
    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions> cost_appro;
    for (auto it = list_.begin(); it != list_.end(); ++it)
    {
      cost_appro += it->getQuadraticApproximation(time, state, input, timeTrajectory, stateTrajectoy, inputTrajectory);
    }
    return cost_appro;
  }

  // add cost to list end
  void add(const StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength> &cost)
  {
    list_.insert(list_.end(), cost);
  }

private:
  IntrusiveList<StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength>> list_;
};
