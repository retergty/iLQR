/**
 * @file Multiplier.hpp
 * @brief 拉格朗日乘子：单约束的惩罚与拉格朗日值，及按时刻的乘子集合与插值。
 */
#pragma once
#include <array>
#include <type_traits>

#include "IntrusiveList.hpp"
#include "LinearInterpolation.hpp"
#include "Types.hpp"

/**
 * @brief 单约束的乘子：惩罚项与拉格朗日项。
 * @tparam Scalar 标量类型。
 */
template <typename Scalar>
struct Multiplier {
  /** @brief 默认构造，penalty 与 lagrangian 为 0。 */
  Multiplier() : Multiplier(0, 0) {}
  /** @brief 用给定惩罚与拉格朗日值构造。 */
  Multiplier(const Scalar penaltyArg, const Scalar lagrangianArg)
      : penalty(penaltyArg), lagrangian(lagrangianArg) {}

  /** @brief 惩罚值。 */
  Scalar penalty;
  /** @brief 拉格朗日值。 */
  Scalar lagrangian;
};

/**
 * @brief 某时刻所有约束的乘子集合：状态等式/不等式、状态-输入等式/不等式。
 * @tparam Scalar 标量类型。
 * @tparam StateEqConstrains 状态等式约束数。
 * @tparam StateIneqConstrains 状态不等式约束数。
 * @tparam StateInputEqConstrains 状态-输入等式约束数。
 * @tparam StateInputIneqConstrains 状态-输入不等式约束数。
 */
template <typename Scalar, int StateEqConstrains, int StateIneqConstrains,
          int StateInputEqConstrains, int StateInputIneqConstrains>
struct MultiplierCollection {
  /** @brief 状态等式约束乘子数组。 */
  std::array<Multiplier<Scalar>, StateEqConstrains> stateEq;
  /** @brief 状态不等式约束乘子数组。 */
  std::array<Multiplier<Scalar>, StateIneqConstrains> stateIneq;
  /** @brief 状态-输入等式约束乘子数组。 */
  std::array<Multiplier<Scalar>, StateInputEqConstrains> stateInputEq;
  /** @brief 状态-输入不等式约束乘子数组。 */
  std::array<Multiplier<Scalar>, StateInputIneqConstrains> stateInputIneq;

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
template <typename Scalar, int StateEqConstrains, int StateIneqConstrains,
          int StateInputEqConstrains, int StateInputIneqConstrains,
          size_t ArrayLen>
Multiplier<Scalar> interpolate(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<Multiplier<Scalar>, ArrayLen>& dataArray) {
  const Scalar penalty =
      interpolate(indexAlpha, dataArray,
                  [](const std::array<Multiplier<Scalar>, ArrayLen>& array,
                     size_t t) -> const Scalar& { return array[t].penalty; });

  const Scalar lagrangian = interpolate(
      indexAlpha, dataArray,
      [](const std::array<Multiplier<Scalar>, ArrayLen>& array,
         size_t t) -> const Scalar& { return array[t].lagrangian; });

  return {penalty, lagrangian};
}

/**
 * @brief 对乘子集合轨迹做线性插值。
 * @param [in] indexAlpha 索引与插值系数对。
 * @param [in] dataArray 乘子集合轨迹。
 * @return 插值得到的乘子集合。
 */
template <typename Scalar, int StateEqConstrains, int StateIneqConstrains,
          int StateInputEqConstrains, int StateInputIneqConstrains,
          size_t ArrayLen>
MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains,
                     StateInputEqConstrains, StateInputIneqConstrains>
interpolate(
    const std::pair<int, Scalar>& indexAlpha,
    const std::array<
        MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains,
                             StateInputEqConstrains, StateInputIneqConstrains>,
        ArrayLen>& dataArray) {
  // number of terms
  MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains,
                       StateInputEqConstrains, StateInputIneqConstrains>
      out;

  // state equality
  for (size_t i = 0; i < StateEqConstrains; i++) {
    Scalar penalty = interpolate(
        indexAlpha, dataArray,
        [i](const std::array<
                MultiplierCollection<
                    Scalar, StateEqConstrains, StateIneqConstrains,
                    StateInputEqConstrains, StateInputIneqConstrains>,
                ArrayLen>& array,
            size_t t) -> const Scalar& { return array[t].stateEq[i].penalty; });
    Scalar lagrangian = interpolate(
        indexAlpha, dataArray,
        [i](const std::array<
                MultiplierCollection<
                    Scalar, StateEqConstrains, StateIneqConstrains,
                    StateInputEqConstrains, StateInputIneqConstrains>,
                ArrayLen>& array,
            size_t t) -> const Scalar& {
          return array[t].stateEq[i].lagrangian;
        });
    out.stateEq[i] = {penalty, lagrangian};
  }  // end of i loop

  // state inequality
  for (size_t i = 0; i < StateIneqConstrains; i++) {
    Scalar penalty = interpolate(
        indexAlpha, dataArray,
        [i](const std::array<
                MultiplierCollection<
                    Scalar, StateEqConstrains, StateIneqConstrains,
                    StateInputEqConstrains, StateInputIneqConstrains>,
                ArrayLen>& array,
            size_t t) -> const Scalar& {
          return array[t].stateIneq[i].penalty;
        });
    Scalar lagrangian = interpolate(
        indexAlpha, dataArray,
        [i](const std::array<
                MultiplierCollection<
                    Scalar, StateEqConstrains, StateIneqConstrains,
                    StateInputEqConstrains, StateInputIneqConstrains>,
                ArrayLen>& array,
            size_t t) -> const Scalar& {
          return array[t].stateIneq[i].lagrangian;
        });
    out.stateIneq[i] = {penalty, lagrangian};
  }  // end of i loop

  // state-input equality
  for (size_t i = 0; i < StateInputEqConstrains; i++) {
    Scalar penalty = interpolate(
        indexAlpha, dataArray,
        [i](const std::array<
                MultiplierCollection<
                    Scalar, StateEqConstrains, StateIneqConstrains,
                    StateInputEqConstrains, StateInputIneqConstrains>,
                ArrayLen>& array,
            size_t t) -> const Scalar& {
          return array[t].stateInputEq[i].penalty;
        });
    Scalar lagrangian = interpolate(
        indexAlpha, dataArray,
        [i](const std::array<
                MultiplierCollection<
                    Scalar, StateEqConstrains, StateIneqConstrains,
                    StateInputEqConstrains, StateInputIneqConstrains>,
                ArrayLen>& array,
            size_t t) -> const Scalar& {
          return array[t].stateInputEq[i].lagrangian;
        });
    out.stateInputEq[i] = {penalty, lagrangian};
  }  // end of i loop

  // state-input inequality
  for (size_t i = 0; i < StateInputIneqConstrains; i++) {
    Scalar penalty = interpolate(
        indexAlpha, dataArray,
        [i](const std::array<
                MultiplierCollection<
                    Scalar, StateEqConstrains, StateIneqConstrains,
                    StateInputEqConstrains, StateInputIneqConstrains>,
                ArrayLen>& array,
            size_t t) -> const Scalar {
          return array[t].stateInputIneq[i].penalty;
        });
    Scalar lagrangian = interpolate(
        indexAlpha, dataArray,
        [i](const std::array<
                MultiplierCollection<
                    Scalar, StateEqConstrains, StateIneqConstrains,
                    StateInputEqConstrains, StateInputIneqConstrains>,
                ArrayLen>& array,
            size_t t) -> const Scalar& {
          return array[t].stateInputIneq[i].lagrangian;
        });
    out.stateInputIneq[i] = {penalty, lagrangian};
  }  // end of i loop

  return out;
}

}  // namespace LinearInterpolation