/**
 * @file Multiplier.hpp
 * @brief 拉格朗日乘子：单约束的惩罚与拉格朗日值，及按时刻的乘子集合与插值。
 */
#pragma once
#include <array>
#include <tuple>
#include <utility>

#include "LinearInterpolation.hpp"
#include "Types.hpp"
#include "iLQRDescriptor.hpp"

/**
 * @brief 单个增广拉格朗日项的乘子：惩罚缩放与向量拉格朗日乘子。
 * @tparam Scalar 标量类型。
 * @tparam CDim 该拉格朗日项内部的约束维度。
 */
template <typename Scalar, int CDim>
struct Multiplier {
  /** @brief 默认构造，penalty 与 lagrangian 均置零。 */
  Multiplier() { setZero(); }

  /** @brief 用给定惩罚缩放与拉格朗日向量构造。 */
  Multiplier(const Scalar penaltyArg, const Vector<Scalar, CDim>& lagrangianArg)
      : penalty(penaltyArg), lagrangian(lagrangianArg) {}

  /** @brief 置零。 */
  void setZero() {
    penalty = Scalar(0);
    lagrangian.setZero();
  }

  /** @brief 惩罚缩放。 */
  Scalar penalty;
  /** @brief 拉格朗日乘子向量。 */
  Vector<Scalar, CDim> lagrangian;
};

template <typename Scalar, typename GroupLayout>
struct MultiplierGroup;

template <typename Scalar, typename... Terms>
struct MultiplierGroup<Scalar, ConstraintGroupLayout<Terms...>> {
  /** @brief 每个 constraint term 对应一个固定维度乘子。 */
  std::tuple<Multiplier<Scalar, Terms::CDim>...> terms;

  /** @brief 与另一 MultiplierGroup 交换内容。 */
  void swap(MultiplierGroup& other) { terms.swap(other.terms); }
};

/**
 * @brief 某时刻所有约束的乘子集合：状态等式/不等式、状态-输入等式/不等式。
 * @tparam Scalar 标量类型。
 * @tparam Layout 乘子布局，需提供各类约束的 ConstraintGroupLayout。
 */
template <typename Scalar, typename Layout>
struct MultiplierCollection {
  using StateEqLayout = typename Layout::StateEqLayout;
  using StateIneqLayout = typename Layout::StateIneqLayout;
  using StateInputEqLayout = typename Layout::StateInputEqLayout;
  using StateInputIneqLayout = typename Layout::StateInputIneqLayout;

  static constexpr int StateEq = StateEqLayout::TotalDim;
  static constexpr int StateIneq = StateIneqLayout::TotalDim;
  static constexpr int StateInputEq = StateInputEqLayout::TotalDim;
  static constexpr int StateInputIneq = StateInputIneqLayout::TotalDim;

  static constexpr int StateEqNumTerms = StateEqLayout::NumTerms;
  static constexpr int StateIneqNumTerms = StateIneqLayout::NumTerms;
  static constexpr int StateInputEqNumTerms = StateInputEqLayout::NumTerms;
  static constexpr int StateInputIneqNumTerms = StateInputIneqLayout::NumTerms;

  /** @brief 状态等式约束乘子组。 */
  MultiplierGroup<Scalar, StateEqLayout> stateEq;
  /** @brief 状态不等式约束乘子组。 */
  MultiplierGroup<Scalar, StateIneqLayout> stateIneq;
  /** @brief 状态-输入等式约束乘子组。 */
  MultiplierGroup<Scalar, StateInputEqLayout> stateInputEq;
  /** @brief 状态-输入不等式约束乘子组。 */
  MultiplierGroup<Scalar, StateInputIneqLayout> stateInputIneq;

  /** @brief 与另一 MultiplierCollection 交换内容。 */
  void swap(MultiplierCollection& other) {
    stateEq.swap(other.stateEq);
    stateIneq.swap(other.stateIneq);
    stateInputEq.swap(other.stateInputEq);
    stateInputIneq.swap(other.stateInputIneq);
  }
};

namespace LinearInterpolation {

/**
 * @brief 对乘子轨迹做线性插值。
 * @param [in] indexAlpha 索引与插值系数 (alpha) 对。
 * @param [in] dataArray 乘子轨迹。
 * @return 插值得到的乘子。
 */
template <typename Scalar, int CDim, size_t ArrayLen>
Multiplier<Scalar, CDim> interpolate(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<Multiplier<Scalar, CDim>, ArrayLen>& dataArray) {
  const Scalar penalty = interpolate(
      indexAlpha, dataArray,
      [](const std::array<Multiplier<Scalar, CDim>, ArrayLen>& array,
         size_t t) -> const Scalar& { return array[t].penalty; });

  const Vector<Scalar, CDim> lagrangian = interpolate(
      indexAlpha, dataArray,
      [](const std::array<Multiplier<Scalar, CDim>, ArrayLen>& array, size_t t)
          -> const Vector<Scalar, CDim>& { return array[t].lagrangian; });

  return {penalty, lagrangian};
}

template <size_t I, typename Scalar, typename GroupLayout, size_t ArrayLen>
auto interpolateGroupTerm(const std::pair<int, Scalar>& indexAlpha,
                          const std::array<MultiplierGroup<Scalar, GroupLayout>,
                                           ArrayLen>& dataArray) {
  using Term = typename GroupLayout::template Term<I>;
  constexpr int CDim = Term::CDim;

  const Scalar penalty = interpolate(
      indexAlpha, dataArray,
      [](const std::array<MultiplierGroup<Scalar, GroupLayout>, ArrayLen>&
             array,
         size_t t) -> const Scalar& {
        return std::get<I>(array[t].terms).penalty;
      });

  const Vector<Scalar, CDim> lagrangian = interpolate(
      indexAlpha, dataArray,
      [](const std::array<MultiplierGroup<Scalar, GroupLayout>, ArrayLen>&
             array,
         size_t t) -> const Vector<Scalar, CDim>& {
        return std::get<I>(array[t].terms).lagrangian;
      });

  return Multiplier<Scalar, CDim>{penalty, lagrangian};
}

template <typename Scalar, typename GroupLayout, size_t ArrayLen, size_t... Is>
MultiplierGroup<Scalar, GroupLayout> interpolateGroupImpl(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<MultiplierGroup<Scalar, GroupLayout>, ArrayLen>& dataArray,
    std::index_sequence<Is...>) {
  MultiplierGroup<Scalar, GroupLayout> out;
  out.terms =
      std::make_tuple(interpolateGroupTerm<Is>(indexAlpha, dataArray)...);
  return out;
}

template <typename Scalar, typename... Terms, size_t ArrayLen>
MultiplierGroup<Scalar, ConstraintGroupLayout<Terms...>> interpolate(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<MultiplierGroup<Scalar, ConstraintGroupLayout<Terms...>>,
                     ArrayLen>& dataArray) {
  using GroupLayout = ConstraintGroupLayout<Terms...>;
  return interpolateGroupImpl(
      indexAlpha, dataArray, std::make_index_sequence<GroupLayout::NumTerms>{});
}

template <size_t I, typename Scalar, typename Layout, size_t ArrayLen,
          typename GroupLayout>
auto interpolateCollectionGroupTerm(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<MultiplierCollection<Scalar, Layout>, ArrayLen>& dataArray,
    const MultiplierGroup<Scalar, GroupLayout>
        MultiplierCollection<Scalar, Layout>::* member) {
  using Term = typename GroupLayout::template Term<I>;
  constexpr int CDim = Term::CDim;

  const Scalar penalty = interpolate(
      indexAlpha, dataArray,
      [member](const std::array<MultiplierCollection<Scalar, Layout>, ArrayLen>&
                   array,
               size_t t) -> const Scalar& {
        return std::get<I>((array[t].*member).terms).penalty;
      });

  const Vector<Scalar, CDim> lagrangian = interpolate(
      indexAlpha, dataArray,
      [member](const std::array<MultiplierCollection<Scalar, Layout>, ArrayLen>&
                   array,
               size_t t) -> const Vector<Scalar, CDim>& {
        return std::get<I>((array[t].*member).terms).lagrangian;
      });

  return Multiplier<Scalar, CDim>{penalty, lagrangian};
}

template <typename Scalar, typename Layout, size_t ArrayLen,
          typename GroupLayout, size_t... Is>
MultiplierGroup<Scalar, GroupLayout> interpolateCollectionGroupImpl(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<MultiplierCollection<Scalar, Layout>, ArrayLen>& dataArray,
    const MultiplierGroup<Scalar, GroupLayout>
        MultiplierCollection<Scalar, Layout>::* member,
    std::index_sequence<Is...>) {
  MultiplierGroup<Scalar, GroupLayout> out;
  out.terms = std::make_tuple(
      interpolateCollectionGroupTerm<Is>(indexAlpha, dataArray, member)...);
  return out;
}

template <typename Scalar, typename Layout, size_t ArrayLen,
          typename GroupLayout>
MultiplierGroup<Scalar, GroupLayout> interpolateCollectionGroup(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<MultiplierCollection<Scalar, Layout>, ArrayLen>& dataArray,
    const MultiplierGroup<Scalar, GroupLayout>
        MultiplierCollection<Scalar, Layout>::* member) {
  return interpolateCollectionGroupImpl(
      indexAlpha, dataArray, member,
      std::make_index_sequence<GroupLayout::NumTerms>{});
}

/**
 * @brief 对乘子集合轨迹做线性插值。
 * @param [in] indexAlpha 索引与插值系数对。
 * @param [in] dataArray 乘子集合轨迹。
 * @return 插值得到的乘子集合。
 */
template <typename Scalar, typename Layout, size_t ArrayLen>
MultiplierCollection<Scalar, Layout> interpolate(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<MultiplierCollection<Scalar, Layout>, ArrayLen>&
        dataArray) {
  using Collection_t = MultiplierCollection<Scalar, Layout>;

  Collection_t out;

  out.stateEq =
      interpolateCollectionGroup(indexAlpha, dataArray, &Collection_t::stateEq);
  out.stateIneq = interpolateCollectionGroup(indexAlpha, dataArray,
                                             &Collection_t::stateIneq);
  out.stateInputEq = interpolateCollectionGroup(indexAlpha, dataArray,
                                                &Collection_t::stateInputEq);
  out.stateInputIneq = interpolateCollectionGroup(
      indexAlpha, dataArray, &Collection_t::stateInputIneq);

  return out;
}

}  // namespace LinearInterpolation