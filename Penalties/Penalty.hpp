#pragma once

#include <tuple>

#include "Approximation/LinearApproximation.hpp"
#include "Approximation/QuadraticApproximation.hpp"
#include "Penalties/AugmentedPenaltyBase.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"

/**
 * @file Penalty.hpp
 * @brief 约束惩罚封装：基于 AugmentedPenaltyBase
 * 计算惩罚值、二次近似及乘子更新。
 *
 * 对约束 h(x,u)，惩罚为 p(t, h, l)；本类用链式法则计算约束-惩罚的二阶近似，
 * 并委托底层惩罚类更新拉格朗日乘子。
 */
/**
 * @brief 单约束惩罚封装：取值、二次近似与乘子初始化/更新，均委托 penalty_ptr_。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 */
template <typename Scalar, int XDim, int UDim>
class Penalty final {
 public:
  /** @brief 用增广惩罚基类指针构造。 */
  Penalty(AugmentedPenaltyBase<Scalar>* penaltyPtr)
      : penalty_ptr_(penaltyPtr) {};

  /** @brief 析构函数。 */
  ~Penalty() = default;

  /** @brief 禁止拷贝。 */
  Penalty(const Penalty& other) = delete;

  /**
   * @brief 获取惩罚代价值 p(t, h, l)。
   * @param [in] t 评估时间。
   * @param [in] h 约束值。
   * @param [in] l 拉格朗日乘子。
   * @return 惩罚值。
   */
  Scalar getValue(const Scalar t, const Scalar h, const Scalar l) const {
    return penalty_ptr_->getValue(t, l, h);
  }

  /**
   * @brief 由约束的线性近似经链式法则得到惩罚的二次近似。
   * @param [in] t 评估时间。
   * @param [in] h 约束的线性近似。
   * @param [in] l 拉格朗日乘子。
   * @return 惩罚的二次近似。
   */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      const Scalar t,
      const ScalarFunctionLinearApproximation<Scalar, XDim, UDim>& h,
      const Scalar l) const {
    Scalar penaltyValue = Scalar(0.0);
    Scalar penaltyDerivative, penaltySecondDerivative;
    std::tie(penaltyValue, penaltyDerivative, penaltySecondDerivative) =
        getPenaltyValue1stDev2ndDev(t, h.f, l);
    const Vector<Scalar, XDim> penaltySecondDev_dhdx =
        penaltySecondDerivative * h.dfdx;

    // 确保仅状态情形下 dfdux 的尺寸正确。
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
        penaltyApproximation;

    penaltyApproximation.f = penaltyValue;
    penaltyApproximation.dfdx = h.dfdx * penaltyDerivative;
    penaltyApproximation.dfdxx = h.dfdx * penaltySecondDev_dhdx.transpose();

    if constexpr (UDim > 0) {
      penaltyApproximation.dfdu = h.dfdu * penaltyDerivative;
      penaltyApproximation.dfdux = h.dfdu * penaltySecondDev_dhdx.transpose();
      penaltyApproximation.dfduu =
          h.dfdu * penaltySecondDerivative * h.dfdu.transpose();
    }

    return penaltyApproximation;
  }

  /**
   * @brief 由约束的二次近似经链式法则得到惩罚的二次近似。
   * @param [in] t 评估时间。
   * @param [in] h 约束的二次近似。
   * @param [in] l 拉格朗日乘子。
   * @return 惩罚的二次近似。
   */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      Scalar t,
      const ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& h,
      const Scalar l) const {
    Scalar penaltyValue = Scalar(0.0);
    Scalar penaltyDerivative, penaltySecondDerivative;
    std::tie(penaltyValue, penaltyDerivative, penaltySecondDerivative) =
        getPenaltyValue1stDev2ndDev(t, h.f, l);
    const Vector<Scalar, XDim> penaltySecondDev_dhdx =
        penaltySecondDerivative * h.dfdx;

    // 确保仅状态情形下 dfdux 的尺寸正确。
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
        penaltyApproximation;

    penaltyApproximation.f = penaltyValue;
    penaltyApproximation.dfdx = h.dfdx * penaltyDerivative;
    penaltyApproximation.dfdxx = h.dfdx * penaltySecondDev_dhdx.transpose();

    penaltyApproximation.dfdxx += penaltyDerivative * h.dfdxx;

    if constexpr (UDim > 0) {
      penaltyApproximation.dfdu = h.dfdu * penaltyDerivative;
      penaltyApproximation.dfdux = h.dfdu * penaltySecondDev_dhdx.transpose();
      penaltyApproximation.dfduu =
          h.dfdu * penaltySecondDerivative * h.dfdu.transpose();

      penaltyApproximation.dfduu += penaltyDerivative * h.dfduu;
      penaltyApproximation.dfdux += penaltyDerivative * h.dfdux;
    }

    return penaltyApproximation;
  }

  /**
   * @brief 根据约束值 h 与当前乘子 l 更新拉格朗日乘子。
   * @param [in] t 时间戳。
   * @param [in] h 约束值。
   * @param [in] l 当前拉格朗日乘子。
   * @return 更新后的乘子。
   */
  Scalar updateMultipliers(Scalar t, const Scalar h, const Scalar l) const {
    return penalty_ptr_->updateMultiplier(t, l, h);
  }

  /** @brief 初始化拉格朗日乘子。 @return 初始乘子。 */
  Scalar initializeMultipliers() const {
    return penalty_ptr_->initializeMultiplier();
  }

 private:
  std::tuple<Scalar, Scalar, Scalar> getPenaltyValue1stDev2ndDev(
      const Scalar t, const Scalar h, const Scalar l) const {
    Scalar penaltyValue = penalty_ptr_->getValue(t, l, h);
    Scalar penaltyDerivative = penalty_ptr_->getDerivative(t, l, h);
    Scalar penaltySecondDerivative = penalty_ptr_->getSecondDerivative(t, l, h);

    return {penaltyValue, penaltyDerivative, penaltySecondDerivative};
  }

  AugmentedPenaltyBase<Scalar>* penalty_ptr_;
};
