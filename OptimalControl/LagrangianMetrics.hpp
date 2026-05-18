/**
 * @file LagrangianMetrics.hpp
 * @brief 拉格朗日指标：单约束项的惩罚值与约束值。
 */
#pragma once
#include "Types.hpp"

/**
 * @brief 单约束项的拉格朗日指标：惩罚值与约束值（标量）。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct LagrangianMetrics {
  /** @brief 默认构造，penalty 与 constraint 为 0。 */
  LagrangianMetrics() : LagrangianMetrics(0, 0) {}
  /** @brief 用给定惩罚与约束值构造。 */
  LagrangianMetrics(Scalar penaltyArg, Scalar constraintArg)
      : penalty(penaltyArg), constraint(constraintArg) {}

  /** @brief 惩罚值。 */
  Scalar penalty;
  /** @brief 约束值。 */
  Scalar constraint;
};