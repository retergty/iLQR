/**
 * @file AugmentedPenaltyBase.hpp
 * @brief
 * 增广惩罚基类接口：约束违反的惩罚项，依赖时间、拉格朗日乘子与约束值，提供取值、一阶/二阶导数及乘子更新。
 */
#pragma once

#include "LinearApproximation.hpp"
#include "QuadraticApproximation.hpp"
#include "Types.hpp"

/**
 * @brief
 * 惩罚函数接口：在代价中加入惩罚项以处理约束违反，假设惩罚为凸；一般为时间、乘子与约束值的函数。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
class AugmentedPenaltyBase {
 public:
  /** @brief 默认构造。 */
  AugmentedPenaltyBase() = default;

  /** 默认析构函数 */
  virtual ~AugmentedPenaltyBase() = default;

  /**
   * 计算某个约束值处的惩罚值。
   *
   * @param [in] t: 约束评估时间。
   * @param [in] l: 拉格朗日乘子。
   * @param [in] h: 约束值。
   * @return 惩罚代价。
   */
  virtual Scalar getValue(const Scalar t, const Scalar l,
                          const Scalar h) const = 0;

  /**
   * 计算某个约束值处的惩罚一阶导数。
   *
   * @param [in] t: 约束评估时间。
   * @param [in] l: 拉格朗日乘子。
   * @param [in] h: 约束值。
   * @return 惩罚对约束值的一阶导数。
   */
  virtual Scalar getDerivative(const Scalar t, const Scalar l,
                               const Scalar h) const = 0;

  /**
   * 计算某个约束值处的惩罚二阶导数。
   *
   * @param [in] t: 约束评估时间。
   * @param [in] l: 拉格朗日乘子。
   * @param [in] h: 约束值。
   * @return 惩罚对约束值的二阶导数。
   */
  virtual Scalar getSecondDerivative(const Scalar t, const Scalar l,
                                     const Scalar h) const = 0;

  /**
   * 更新拉格朗日乘子。
   *
   * @param [in] t: 时间戳。
   * @param [in] l: 拉格朗日乘子。
   * @param [in] h: 约束值。
   * @return 更新后的拉格朗日乘子。
   */
  virtual Scalar updateMultiplier(const Scalar t, const Scalar l,
                                  const Scalar h) const = 0;

  /**
   * 初始化拉格朗日乘子。
   *
   * @return 拉格朗日乘子的初始值。
   */
  virtual Scalar initializeMultiplier() const = 0;

  AugmentedPenaltyBase(const AugmentedPenaltyBase& other) = default;
};
