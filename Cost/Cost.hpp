/**
 * @file Cost.hpp
 * @brief 代价项接口：仅状态代价与状态-输入代价，提供取值与二次近似。
 */
#pragma once
#include "Types.hpp"
#include "Reference.hpp"
#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "IntrusiveList.hpp"
#include "LinearInterpolation.hpp"

/**
 * @brief 仅状态代价项基类：按时间与状态计算代价值及二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDimisions, int ArrayLength>
class StateCost : public IntrusiveListNode<StateCost<Scalar, XDimisions, ArrayLength>>
{
public:
  /** @brief 构造，cost_number 为代价项唯一标识。 */
  StateCost(int cost_number) : number(cost_number) {};
  virtual ~StateCost() = default;

  /** @brief 判断该时刻代价项是否激活。 */
  virtual bool isActive(Scalar time) const { 
      (void)time;
      return true; }

  /** @brief 获取代价值。 */
  virtual Scalar getValue(Scalar time, const Vector<Scalar, XDimisions>& state,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const = 0;
  virtual Scalar getValue(
      int time_index, const Vector<Scalar, XDimisions>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const = 0;

  /** @brief 获取代价的二次近似（仅状态）。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
    getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions>& state,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const = 0;

  /** @brief 获取代价的二次近似（仅状态）。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDimisions, 0>
      getQuadraticApproximation(
          int time_index, const Vector<Scalar, XDimisions>& state,
          const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy) const = 0;

  /** @brief 代价项唯一标识号。 */
  int number;

protected:
  StateCost(const StateCost& rhs) = default;
};

/**
 * @brief 状态-输入代价项基类：按时间、状态与输入计算代价值及二次近似。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam UDimisions 输入维度。
 * @tparam ArrayLength 轨迹长度。
 */
template <typename Scalar, int XDimisions, int UDimisions, int ArrayLength>
class StateInputCost : public IntrusiveListNode<StateInputCost<Scalar, XDimisions, UDimisions, ArrayLength>>
{
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
  virtual Scalar getValue(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const = 0;
  virtual Scalar getValue(int time_index, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
    const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const = 0;

  /** @brief 获取代价的二次近似（状态-输入）。 */
  virtual ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>
    getQuadraticApproximation(Scalar time, const Vector<Scalar, XDimisions>& state, const Vector<Scalar, UDimisions>& input,
      const std::array<Scalar, ArrayLength>& timeTrajectory, const std::array<Vector<Scalar, XDimisions>, ArrayLength>& stateTrajectoy, const std::array<Vector<Scalar, UDimisions>, ArrayLength>& inputTrajectory) const = 0;

  // identify state input cost, must be unique
  int number;

protected:
  StateInputCost(const StateInputCost& rhs) = default;
};