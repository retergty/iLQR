/**
 * @file CostCollection.hpp
 * @brief 代价项集合：汇总多个仅状态或状态-输入代价项，提供总代价与总二次近似。
 */
#pragma once

#include "Cost/Cost.hpp"
#include "IntrusiveList/IntrusiveList.hpp"

/**
 * @brief 仅状态代价项集合：对多个 StateCost 求和得到总代价与总二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDim, int UDim, int ArrayLength>
class StateCostCollection {
 public:
  StateCostCollection() = default;
  virtual ~StateCostCollection() = default;

  /** @brief 获取仅状态总代价值。 */
  Scalar getValue(
      Scalar time, const Vector<Scalar, XDim>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectories,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoies) {
    Scalar cost = 0;
    for (auto it = list_.begin(); it != list_.end(); ++it) {
      if (it->isActive(time)) {
        cost += it->getValue(time, state, timeTrajectories, stateTrajectoies);
      }
    }
    return cost;
  }

  /**
   * @brief 累加仅状态总代价的二次近似。
   * @note 本接口不会清零 addAppro；调用方需自行管理初始化。
   */
  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectories,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoies,
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& addAppro) {
    for (auto it = list_.begin(); it != list_.end(); ++it) {
      if (it->isActive(time)) {
        it->addQuadraticApproximation(time, state, timeTrajectories,
                                      stateTrajectoies, addAppro);
      }
    }
  }

  /** @brief 在链表末尾添加代价项。 */
  void add(StateCost<Scalar, XDim, UDim, ArrayLength>& cost) {
    list_.insert(list_.end(), cost);
  }

 private:
  IntrusiveList<StateCost<Scalar, XDim, UDim, ArrayLength>> list_;
};

/**
 * @brief 状态-输入代价项集合：对多个 StateInputCost
 * 求和得到总代价与总二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDim, int UDim, int ArrayLength>
class StateInputCostCollection {
 public:
  StateInputCostCollection() = default;
  ~StateInputCostCollection() = default;

  /** 获取状态-输入代价值。 */
  Scalar getValue(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory) {
    Scalar cost = 0;
    for (auto it = list_.begin(); it != list_.end(); ++it) {
      if (it->isActive(time)) {
        cost += it->getValue(time, state, input, timeTrajectory, stateTrajectoy,
                             inputTrajectory);
      }
    }
    return cost;
  }

  /**
   * @brief 累加状态-输入代价二次近似。
   * @note 本接口不会清零 addAppro；调用方需自行管理初始化。
   */
  void addQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory,
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& addAppro) {
    for (auto it = list_.begin(); it != list_.end(); ++it) {
      if (it->isActive(time)) {
        it->addQuadraticApproximation(time, state, input, timeTrajectory,
                                      stateTrajectoy, inputTrajectory,
                                      addAppro);
      }
    }
  }

  // 将代价加入列表尾部。
  void add(StateInputCost<Scalar, XDim, UDim, ArrayLength>& cost) {
    list_.insert(list_.end(), cost);
  }

 private:
  IntrusiveList<StateInputCost<Scalar, XDim, UDim, ArrayLength>> list_;
};
