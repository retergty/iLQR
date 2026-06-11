/**
 * @file StateInputAugmentedLagrangianCollection.hpp
 * @brief
 * 状态-输入增广拉格朗日集合：汇总多个状态-输入约束的惩罚项，提供总取值与总二次近似。
 */
#pragma once

#include <cassert>
#include <tuple>
#include <utility>

#include "AugmentedLagrangian/StateInputAugmentedLagrangianInterface.hpp"
#include "ModelData/Metrics.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"

/**
 * @brief 状态-输入增广拉格朗日惩罚项集合：对多个 StateInputAugmentedLagrangian
 * term 求和。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 输入维度。
 * @tparam GroupLayout 约束 term 分组布局。
 */
template <typename Scalar, int XDim, int UDim, typename GroupLayout>
class StateInputAugmentedLagrangianCollection;

template <typename Scalar, int XDim, int UDim, typename... Terms>
class StateInputAugmentedLagrangianCollection<Scalar, XDim, UDim,
                                              ConstraintGroupLayout<Terms...>> {
 public:
  using Layout = ConstraintGroupLayout<Terms...>;

  StateInputAugmentedLagrangianCollection() = default;

  template <std::size_t I>
  using TermLayout = typename Layout::template Term<I>;

  template <std::size_t I>
  using TermPtr =
      const StateInputAugmentedLagrangianInterface<Scalar, XDim, UDim,
                                                   TermLayout<I>::CDim>*;

  /** @brief 设置第 I 个增广拉格朗日项。 */
  template <std::size_t I>
  void set(TermPtr<I> term) {
    assert(term != nullptr);
    std::get<I>(terms_) = term;
  }

  /** @brief 获取各 term 的状态-输入约束与惩罚值。 */
  LagrangianMetricsGroup<Scalar, Layout> getValue(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const MultiplierGroup<Scalar, Layout>& termsMultiplier) const {
    return getValueImpl(time, state, input, termsMultiplier,
                        std::make_index_sequence<Layout::NumTerms>{});
  }

  /** @brief 计算并累加所有 state-input Lagrangian 惩罚项的二次近似。 */
  void addQuadraticApproximation(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const MultiplierGroup<Scalar, Layout>& termsMultiplier,
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& addAppro)
      const {
    addQuadraticApproximationImpl(
        time, state, input, termsMultiplier, addAppro,
        std::make_index_sequence<Layout::NumTerms>{});
  }

  /** @brief 获取所有 state-input Lagrangian 惩罚项的二次近似之和。 */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const MultiplierGroup<Scalar, Layout>& termsMultiplier) const {
    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> penalty;
    penalty.setZero();
    addQuadraticApproximation(time, state, input, termsMultiplier, penalty);
    return penalty;
  }

  /** 更新每个激活项的拉格朗日/惩罚乘子和惩罚值。
   * 项。 */
  void updateLagrangian(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      LagrangianMetricsGroup<Scalar, Layout>& termsMetrics,
      MultiplierGroup<Scalar, Layout>& termsMultiplier) const {
    updateLagrangianImpl(time, state, input, termsMetrics, termsMultiplier,
                         std::make_index_sequence<Layout::NumTerms>{});
  }

  /** 初始化每个激活项的拉格朗日/惩罚乘子。 */
  void initializeLagrangian(
      const Scalar time,
      MultiplierGroup<Scalar, Layout>& termsMultiplier) const {
    initializeLagrangianImpl(time, termsMultiplier,
                             std::make_index_sequence<Layout::NumTerms>{});
  }

 private:
  template <std::size_t I>
  TermPtr<I> getTerm() const {
    const auto* term = std::get<I>(terms_);
    assert(term != nullptr);
    return term;
  }

  template <std::size_t... Is>
  LagrangianMetricsGroup<Scalar, Layout> getValueImpl(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const MultiplierGroup<Scalar, Layout>& termsMultiplier,
      std::index_sequence<Is...>) const {
    LagrangianMetricsGroup<Scalar, Layout> result;
    result.terms = std::make_tuple(getTerm<Is>()->getValue(
        time, state, input, std::get<Is>(termsMultiplier.terms))...);
    return result;
  }

  template <std::size_t... Is>
  void addQuadraticApproximationImpl(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const MultiplierGroup<Scalar, Layout>& termsMultiplier,
      ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>& addAppro,
      std::index_sequence<Is...>) const {
    ((getTerm<Is>()->addQuadraticApproximation(
          time, state, input, std::get<Is>(termsMultiplier.terms), addAppro)),
     ...);
  }

  template <std::size_t... Is>
  void updateLagrangianImpl(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      LagrangianMetricsGroup<Scalar, Layout>& termsMetrics,
      MultiplierGroup<Scalar, Layout>& termsMultiplier,
      std::index_sequence<Is...>) const {
    (void)time;
    (void)state;
    (void)input;
    ((std::tie(std::get<Is>(termsMultiplier.terms),
               std::get<Is>(termsMetrics.terms).penalty) =
          getTerm<Is>()->updateLagrangian(
              time, state, input, std::get<Is>(termsMetrics.terms).constraint,
              std::get<Is>(termsMultiplier.terms))),
     ...);
  }

  template <std::size_t... Is>
  void initializeLagrangianImpl(
      const Scalar time, MultiplierGroup<Scalar, Layout>& termsMultiplier,
      std::index_sequence<Is...>) const {
    (void)time;
    ((std::get<Is>(termsMultiplier.terms) =
          getTerm<Is>()->initializeLagrangian(time)),
     ...);
  }

  std::tuple<const StateInputAugmentedLagrangianInterface<Scalar, XDim, UDim,
                                                          Terms::CDim>*...>
      terms_;
};