/**
 * @file CostCollection.hpp
 * @brief 代价项集合：汇总多个仅状态或状态-输入代价项，提供总代价与总二次近似。
 */
#pragma once
#include "Cost.hpp"
#include <memory>
#include "IntrusiveList.hpp"

/**
 * @brief 仅状态代价项集合：对多个 StateCost 求和得到总代价与总二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDimisions, int ArrayLength>
class StateCostCollection
{
public:
  StateCostCollection() = default;
  virtual ~StateCostCollection() = default;

  /** @brief 获取仅状态总代价值。 */
  Scalar getValue(Scalar time, const Vector<Scalar, XDimisions> &state, const std::array<Scalar, ArrayLength> &timeTrajectories, const std::array<Vector<Scalar, XDimisions>, ArrayLength> &stateTrajectoies) const
  {
    Scalar cost = 0;
    for (auto it = list_.begin(); it != list_.end(); ++it)
    {
      cost += it->getValue(time, state, timeTrajectories, stateTrajectoies);
    }
    return cost;
  }

  /** @brief 获取仅状态总代价的二次近似。 */
  ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
  getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions> &state, const std::array<Scalar, ArrayLength> &timeTrajectories, const std::array<Vector<Scalar, XDimisions>, ArrayLength> &stateTrajectoies) const
  {
    ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0> cost_appro;
    cost_appro.setZero();
    for (auto it = list_.begin(); it != list_.end(); ++it)
    {
      cost_appro += it->getQuadraticApproximation(time, state, timeTrajectories, stateTrajectoies);
    }
    return cost_appro;
  }

  /** @brief 在链表末尾添加代价项。 */
  void add(StateCost<Scalar, XDimisions, ArrayLength> &cost)
  {
    list_.insert(list_.end(), cost);
  }

private:
  IntrusiveList<StateCost<Scalar, XDimisions, ArrayLength>> list_;
};

/**
 * @brief 状态-输入代价项集合：对多个 StateInputCost 求和得到总代价与总二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam UDimisions 输入维度。
 * @tparam ArrayLength 轨迹长度。
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
    cost_appro.setZero();
    for (auto it = list_.begin(); it != list_.end(); ++it)
    {
      cost_appro += it->getQuadraticApproximation(time, state, input, timeTrajectory, stateTrajectoy, inputTrajectory);
    }
    return cost_appro;
  }

  // add cost to list end
  void add(StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength> &cost)
  {
    list_.insert(list_.end(), cost);
  }

private:
  IntrusiveList<StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength>> list_;
};
