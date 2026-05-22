/**
 * @file LagrangianMetrics.hpp
 * @brief 拉格朗日指标：单个增广拉格朗日项的惩罚值与约束违反值。
 */
#pragma once
#include "Types.hpp"

/**
 * @brief 单个增广拉格朗日项的指标：惩罚值与约束违反值。
 *
 * 当前实现仍保存标量约束值；在向量约束重构完成后，`constraint`
 * 将扩展为固定维度向量，用于表示同一个拉格朗日项内的多个约束分量。
 *
 * @tparam Scalar 标量类型。
 */
template <typename Scalar, int CDim>
struct LagrangianMetrics {
  LagrangianMetrics() { setZero(); }
  LagrangianMetrics(Scalar penaltyArg,
                    const Vector<Scalar, CDim>& constraintArg)
      : penalty(penaltyArg), constraint(constraintArg) {}
  void setZero() {
    penalty = Scalar(0);
    constraint.setZero();
  }
  Scalar penalty;
  Vector<Scalar, CDim> constraint;
};