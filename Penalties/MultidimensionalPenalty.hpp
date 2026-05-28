#pragma once

#include <cassert>
#include <tuple>

#include "Approximation/LinearApproximation.hpp"
#include "Approximation/QuadraticApproximation.hpp"
#include "Matrix/Types.hpp"
#include "Penalties/AugmentedPenaltyBase.hpp"

/**
 * @brief 多维约束惩罚封装。
 *
 * 对向量约束 \f$ h(x,u) \in \mathbb{R}^{CDim} \f$，该类逐分量复用当前
 * 标量 `AugmentedPenaltyBase`，并将各分量惩罚求和：
 *
 * \f[
 *   p(t, h, \lambda) = \sum_i p_i(t, h_i, \lambda_i)
 * \f]
 *
 * 该类负责使用链式法则，把向量约束的一阶/二阶近似转换为最终标量
 * cost 的二次近似。
 */
template <typename Scalar, int XDim, int UDim, int CDim>
class MultidimensionalPenalty final {
 public:
  /** @brief 用一个标量惩罚函数构造；该惩罚会应用到每个约束分量。 */
  MultidimensionalPenalty(const AugmentedPenaltyBase<Scalar>* penaltyPtr)
      : penalty_ptr_(penaltyPtr) {
    assert(penalty_ptr_ != nullptr);
  }

  /** @brief 默认析构函数。 */
  ~MultidimensionalPenalty() = default;

  /**
   * @brief 获取向量约束的总惩罚值。
   * @param [in] t 评估时间。
   * @param [in] h 约束值向量。
   * @param [in] l 拉格朗日乘子向量。
   * @return 各约束分量惩罚值之和。
   */
  Scalar getValue(const Scalar t, const Vector<Scalar, CDim>& h,
                  const Vector<Scalar, CDim>& l) const {
    Scalar penalty = 0;
    for (int i = 0; i < CDim; ++i) {
      penalty += penalty_ptr_->getValue(t, l(i), h(i));
    }

    return penalty;
  }

  /**
   * @brief 由向量约束的线性近似经链式法则得到标量惩罚的二次近似。
   * @param [in] t 评估时间。
   * @param [in] h 约束的线性近似。
   * @param [in] l 拉格朗日乘子向量。
   * @return 标量惩罚函数对状态/输入的二次近似。
   */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      const Scalar t,
      const VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim>& h,
      const Vector<Scalar, CDim>& l) const {
    Scalar penaltyValue = 0.0;
    Vector<Scalar, CDim> penaltyDerivative, penaltySecondDerivative;
    std::tie(penaltyValue, penaltyDerivative, penaltySecondDerivative) =
        getPenaltyValue1stDev2ndDev(t, h.f, l);
    const Matrix<Scalar, CDim, XDim> penaltySecondDev_dhdx =
        scaleRows(penaltySecondDerivative, h.dfdx);

    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
        penaltyApproximation;
    penaltyApproximation.setZero();

    penaltyApproximation.f = penaltyValue;
    penaltyApproximation.dfdx = h.dfdx.transpose() * penaltyDerivative;
    penaltyApproximation.dfdxx = h.dfdx.transpose() * penaltySecondDev_dhdx;
    if constexpr (UDim > 0) {
      penaltyApproximation.dfdu = h.dfdu.transpose() * penaltyDerivative;
      penaltyApproximation.dfdux = h.dfdu.transpose() * penaltySecondDev_dhdx;
      penaltyApproximation.dfduu =
          h.dfdu.transpose() * scaleRows(penaltySecondDerivative, h.dfdu);
    }

    return penaltyApproximation;
  }

  /**
   * @brief 由向量约束的二次近似经链式法则得到标量惩罚的二次近似。
   * @param [in] t 评估时间。
   * @param [in] h 约束的二次近似。
   * @param [in] l 拉格朗日乘子向量。
   * @return 标量惩罚函数对状态/输入的二次近似。
   */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      const Scalar t,
      const VectorFunctionQuadraticApproximation<Scalar, CDim, XDim, UDim>& h,
      const Vector<Scalar, CDim>& l) const {
    Scalar penaltyValue = 0.0;
    Vector<Scalar, CDim> penaltyDerivative, penaltySecondDerivative;
    std::tie(penaltyValue, penaltyDerivative, penaltySecondDerivative) =
        getPenaltyValue1stDev2ndDev(t, h.f, l);
    const Matrix<Scalar, CDim, XDim> penaltySecondDev_dhdx =
        scaleRows(penaltySecondDerivative, h.dfdx);

    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
        penaltyApproximation;
    penaltyApproximation.setZero();

    penaltyApproximation.f = penaltyValue;
    penaltyApproximation.dfdx = h.dfdx.transpose() * penaltyDerivative;
    penaltyApproximation.dfdxx = h.dfdx.transpose() * penaltySecondDev_dhdx;
    for (int i = 0; i < CDim; ++i) {
      penaltyApproximation.dfdxx += penaltyDerivative(i) * h.dfdxx[i];
    }

    if constexpr (UDim > 0) {
      penaltyApproximation.dfdu = h.dfdu.transpose() * penaltyDerivative;
      penaltyApproximation.dfdux = h.dfdu.transpose() * penaltySecondDev_dhdx;
      penaltyApproximation.dfduu =
          h.dfdu.transpose() * scaleRows(penaltySecondDerivative, h.dfdu);
      for (int i = 0; i < CDim; ++i) {
        penaltyApproximation.dfduu += penaltyDerivative(i) * h.dfduu[i];
        penaltyApproximation.dfdux += penaltyDerivative(i) * h.dfdux[i];
      }
    }

    return penaltyApproximation;
  }

  /**
   * @brief 更新拉格朗日乘子向量。
   * @param [in] t 时间戳。
   * @param [in] h 约束值向量。
   * @param [in] l 当前拉格朗日乘子向量。
   * @return 更新后的拉格朗日乘子向量。
   */
  Vector<Scalar, CDim> updateMultipliers(const Scalar t,
                                         const Vector<Scalar, CDim>& h,
                                         const Vector<Scalar, CDim>& l) const {
    Vector<Scalar, CDim> updated_l;
    for (int i = 0; i < CDim; ++i) {
      updated_l(i) = penalty_ptr_->updateMultiplier(t, l(i), h(i));
    }

    return updated_l;
  }

  /**
   * @brief 初始化拉格朗日乘子向量。
   * @return 每个约束分量的初始拉格朗日乘子。
   */
  Vector<Scalar, CDim> initializeMultipliers() const {
    Vector<Scalar, CDim> l;
    for (int i = 0; i < CDim; ++i) {
      l(i) = penalty_ptr_->initializeMultiplier();
    }

    return l;
  }

 private:
  template <int Cols>
  Matrix<Scalar, CDim, Cols> scaleRows(
      const Vector<Scalar, CDim>& diagonal,
      const Matrix<Scalar, CDim, Cols>& matrix) const {
    Matrix<Scalar, CDim, Cols> result;
    for (int i = 0; i < CDim; ++i) {
      for (int j = 0; j < Cols; ++j) {
        result(i, j) = diagonal(i) * matrix(i, j);
      }
    }
    return result;
  }

  std::tuple<Scalar, Vector<Scalar, CDim>, Vector<Scalar, CDim>>
  getPenaltyValue1stDev2ndDev(Scalar t, const Vector<Scalar, CDim>& h,
                              const Vector<Scalar, CDim>& l) const {
    Scalar penaltyValue = 0.0;
    Vector<Scalar, CDim> penaltyDerivative;
    Vector<Scalar, CDim> penaltySecondDerivative;
    for (int i = 0; i < CDim; ++i) {
      penaltyValue += penalty_ptr_->getValue(t, l(i), h(i));
      penaltyDerivative(i) = penalty_ptr_->getDerivative(t, l(i), h(i));
      penaltySecondDerivative(i) =
          penalty_ptr_->getSecondDerivative(t, l(i), h(i));
    }

    return {penaltyValue, penaltyDerivative, penaltySecondDerivative};
  }

  const AugmentedPenaltyBase<Scalar>* penalty_ptr_;
};
