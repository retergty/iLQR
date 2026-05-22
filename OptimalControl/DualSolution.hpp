/**
 * @file DualSolution.hpp
 * @brief 对偶解：各时刻的乘子集合及按时间插值查询。
 */
#pragma once
#include <array>

#include "LinearInterpolation.hpp"
#include "Multiplier.hpp"
#include "Numerics.hpp"
#include "iLQRDescriptor.hpp"

/**
 * @brief 对偶解：时间序列、终端乘子与各中间时刻乘子集合。
 * @tparam Scalar 标量类型。
 * @tparam Horizon 预测时域配置，提供 PredictLength。
 * @tparam ConstraintConfig 约束配置。
 */
template <typename Scalar, typename Horizon, typename ConstraintConfig>
struct DualSolution {
  static constexpr std::size_t PredictLength = Horizon::PredictLength;

  using IntermediateMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           IntermediateStageConstraintLayout<ConstraintConfig>>;
  using FinalMultiplierCollection_t =
      MultiplierCollection<Scalar,
                           FinalStageConstraintLayout<ConstraintConfig>>;

  /** @brief 默认构造为空缓存状态。 */
  DualSolution() { clear(); }

  /** @brief 时间序列，与轨迹对齐。 */
  std::array<Scalar, PredictLength + 1> timeTrajectory;

  /** @brief 终端时刻乘子集合。 */
  FinalMultiplierCollection_t final;
  /** @brief 各中间时刻乘子集合。 */
  std::array<IntermediateMultiplierCollection_t, PredictLength> intermediates;

  /** @brief 与另一 DualSolution 交换内容。 */
  void swap(DualSolution& other) {
    timeTrajectory.swap(other.timeTrajectory);
    intermediates.swap(other.intermediates);
    final.swap(other.final);
  }

  /** @brief 判断是否为空：时间序列是否全为 0。 */
  bool empty() const {
    for (size_t i = 0; i < timeTrajectory.size(); ++i) {
      if (!numerics::almost_eq(timeTrajectory[i], static_cast<Scalar>(0))) {
        return false;
      }
    }
    return true;
  }

  /** @brief 清空：时间序列置零。 */
  void clear() {
    for (size_t i = 0; i < timeTrajectory.size(); ++i) {
      timeTrajectory[i] = 0;
    }
  }
};

/**
 * @brief 对偶解的引用视图：仅引用 terminal 与
 * intermediates，不包含时间戳，表示不可修改时间戳。
 */
template <typename Scalar, typename Horizon, typename ConstraintConfig>
struct DualSolutionRef {
  static constexpr std::size_t PredictLength = Horizon::PredictLength;

  using DualSolution_t = DualSolution<Scalar, Horizon, ConstraintConfig>;
  using IntermediateMultiplierCollection_t =
      typename DualSolution_t::IntermediateMultiplierCollection_t;
  using FinalMultiplierCollection_t =
      typename DualSolution_t::FinalMultiplierCollection_t;

  /** @brief 由 DualSolution 构造，引用其 final 与 intermediates。 */
  DualSolutionRef(DualSolution_t& dualSolution)
      : DualSolutionRef(dualSolution.final, dualSolution.intermediates) {}

  /** @brief 直接引用给定的 final 与 intermediates 数组。 */
  DualSolutionRef(FinalMultiplierCollection_t& finalRef,
                  std::array<IntermediateMultiplierCollection_t, PredictLength>&
                      intermediatesRef)
      : final(finalRef), intermediates(intermediatesRef) {}

  /** @brief 终端乘子集合的引用。 */
  FinalMultiplierCollection_t& final;
  /** @brief 中间时刻乘子集合的引用。 */
  std::array<IntermediateMultiplierCollection_t, PredictLength>& intermediates;
};

/**
 * @brief 按给定时间在对偶解时间序列上插值，得到该时刻的中间乘子集合。
 * @param [in] dualSolution 对偶解。
 * @param [in] time 查询时间。
 * @return 该时刻对应的状态/状态-输入、等式/不等式拉格朗日乘子集合（插值结果）。
 */
template <typename Scalar, typename Horizon, typename ConstraintConfig>
inline
    typename DualSolution<Scalar, Horizon,
                          ConstraintConfig>::IntermediateMultiplierCollection_t
    getIntermediateDualSolutionAtTime(
        const DualSolution<Scalar, Horizon, ConstraintConfig>& dualSolution,
        Scalar time) {
  const std::pair<int, Scalar> indexAlpha =
      LinearInterpolation::timeSegment(time, dualSolution.timeTrajectory);
  return LinearInterpolation::interpolate(indexAlpha,
                                          dualSolution.intermediates);
}