/**
 * @file StateAugmentedLagrangianCollection.hpp
 * @brief
 * 仅状态增广拉格朗日集合：汇总多个仅状态约束的惩罚项，提供总取值与总二次近似。
 */
#pragma once

#include <cassert>
#include <tuple>
#include <utility>

#include "AugmentedLagrangian/StateAugmentedLagrangianInterface.hpp"
#include "ModelData/Metrics.hpp"
#include "iLQR/LinearAlgebraTypes.hpp"

/**
 * @brief 仅状态增广拉格朗日惩罚项集合：对多个不同维度的
 * StateAugmentedLagrangian term 求和。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam GroupLayout 约束 term 分组布局。
 */
template <typename Scalar, int XDim, typename GroupLayout>
class StateAugmentedLagrangianCollection;

template <typename Scalar, int XDim, typename... Terms>
class StateAugmentedLagrangianCollection<Scalar, XDim,
                                         ConstraintGroupLayout<Terms...>> {
 public:
  using Layout = ConstraintGroupLayout<Terms...>;

  StateAugmentedLagrangianCollection() = default;

  template <std::size_t I>
  using TermLayout = typename Layout::template Term<I>;

  template <std::size_t I>
  using TermPtr = const StateAugmentedLagrangianInterface<Scalar, XDim,
                                                          TermLayout<I>::CDim>*;

  /** @brief 设置第 I 个增广拉格朗日项。 */
  template <std::size_t I>
  void set(TermPtr<I> term) {
    assert(term != nullptr);
    std::get<I>(terms_) = term;
  }

  /** @brief 获取各 term 的约束与惩罚值。 */
  LagrangianMetricsGroup<Scalar, Layout> getValue(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const MultiplierGroup<Scalar, Layout>& termsMultiplier) const {
    return getValueImpl(time, state, termsMultiplier,
                        std::make_index_sequence<Layout::NumTerms>{});
  }

  /** @brief 获取所有 state Lagrangian 惩罚项的二次近似之和。 */
  ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>
  getQuadraticApproximation(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const MultiplierGroup<Scalar, Layout>& termsMultiplier) const {
    ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> penalty;
    penalty.setZero();
    addQuadraticApproximationImpl(
        time, state, termsMultiplier, penalty,
        std::make_index_sequence<Layout::NumTerms>{});
    return penalty;
  }

  /** 更新每个激活项的拉格朗日/惩罚乘子以及惩罚值。
   * 项。 */
  void updateLagrangian(
      Scalar time, const Vector<Scalar, XDim>& state,
      LagrangianMetricsGroup<Scalar, Layout>& termsMetrics,
      MultiplierGroup<Scalar, Layout>& termsMultiplier) const {
    updateLagrangianImpl(time, state, termsMetrics, termsMultiplier,
                         std::make_index_sequence<Layout::NumTerms>{});
  }

  /** 初始化每个激活项的拉格朗日/惩罚乘子。 */
  void initializeLagrangian(
      Scalar time, MultiplierGroup<Scalar, Layout>& termsMultiplier) const {
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
      const MultiplierGroup<Scalar, Layout>& termsMultiplier,
      std::index_sequence<Is...>) const {
    LagrangianMetricsGroup<Scalar, Layout> result;
    result.terms = std::make_tuple(getTerm<Is>()->getValue(
        time, state, std::get<Is>(termsMultiplier.terms))...);
    return result;
  }

  template <std::size_t... Is>
  void addQuadraticApproximationImpl(
      const Scalar time, const Vector<Scalar, XDim>& state,
      const MultiplierGroup<Scalar, Layout>& termsMultiplier,
      ScalarFunctionQuadraticApproximation<Scalar, XDim, 0>& addAppro,
      std::index_sequence<Is...>) const {
    ((getTerm<Is>()->addQuadraticApproximation(
          time, state, std::get<Is>(termsMultiplier.terms), addAppro)),
     ...);
  }

  template <std::size_t... Is>
  void updateLagrangianImpl(
      const Scalar time, const Vector<Scalar, XDim>& state,
      LagrangianMetricsGroup<Scalar, Layout>& termsMetrics,
      MultiplierGroup<Scalar, Layout>& termsMultiplier,
      std::index_sequence<Is...>) const {
    (void)time;
    (void)state;
    ((std::tie(std::get<Is>(termsMultiplier.terms),
               std::get<Is>(termsMetrics.terms).penalty) =
          getTerm<Is>()->updateLagrangian(
              time, state, std::get<Is>(termsMetrics.terms).constraint,
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

  std::tuple<
      const StateAugmentedLagrangianInterface<Scalar, XDim, Terms::CDim>*...>
      terms_;
};
