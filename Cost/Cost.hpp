/**
 * @file Cost.hpp
 * @brief 代价项接口：仅状态代价与状态-输入代价，提供取值与二次近似。
 */
#pragma once
#include "IntrusiveList.hpp"
#include "QuadraticApproximation.hpp"
#include "Types.hpp"

/**
 * @brief 仅状态代价项基类：按时间与状态计算代价值及二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDim, int ArrayLength>
class StateCost
    : public IntrusiveListNode<StateCost<Scalar, XDim, ArrayLength>> {
 public:
  /** @brief 构造，cost_number 为代价项唯一标识。 */
  StateCost(int cost_number) : number(cost_number) {};
  virtual ~StateCost() = default;

  /** @brief 判断该时刻代价项是否激活。 */
  virtual bool isActive(Scalar time) const {
    (void)time;
    return true;
  }

  /** @brief 获取代价值。 */
  virtual Scalar getValue(Scalar time, const Vector<Scalar, XDim>& state,
                          const std::array<Scalar, ArrayLength>& timeTrajectory,
                          const std::array<Vector<Scalar, XDim>, ArrayLength>&
                              stateTrajectoy) const = 0;
  virtual Scalar getValue(int time_index, const Vector<Scalar, XDim>& state,
                          const std::array<Scalar, ArrayLength>& timeTrajectory,
                          const std::array<Vector<Scalar, XDim>, ArrayLength>&
                              stateTrajectoy) const = 0;

  /** @brief 获取代价的二次近似（仅状态）。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy)
      const = 0;

  /** @brief 获取代价的二次近似（仅状态）。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>
  getQuadraticApproximation(
      int time_index, const Vector<Scalar, XDim>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy)
      const = 0;

  /** @brief 代价项唯一标识号。 */
  int number;

 protected:
  StateCost(const StateCost& rhs) = default;
};

/**
 * @brief 状态-输入代价项基类：按时间、状态与输入计算代价值及二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDim, int UDim, int ArrayLength>
class StateInputCost : public IntrusiveListNode<
                           StateInputCost<Scalar, XDim, UDim, ArrayLength>> {
 public:
  /** @brief 构造，cost_number 为代价项唯一标识。 */
  StateInputCost(int cost_number) : number(cost_number) {};
  virtual ~StateInputCost() = default;

  /** @brief 判断该时刻代价项是否激活。 */
  virtual bool isActive(Scalar time) const {
    (void)time;
    return true;
  }

  /** @brief 获取代价值。 */
  virtual Scalar getValue(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory)
      const = 0;
  virtual Scalar getValue(
      int time_index, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory)
      const = 0;

  /** @brief 获取代价的二次近似（状态-输入）。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory,
      const std::array<Vector<Scalar, XDim>, ArrayLength>& stateTrajectoy,
      const std::array<Vector<Scalar, UDim>, ArrayLength>& inputTrajectory)
      const = 0;

  // 标识状态-输入代价，必须唯一。
  int number;

 protected:
  StateInputCost(const StateInputCost& rhs) = default;
};