#pragma once

#include "Types.hpp"

namespace qp_solver {

/** 时间、状态、输入轨迹。最后一个时间点只有状态，没有
 * 输入。 */
template <typename Scalar, int XDim, int UDim, size_t PredictLength>
struct ContinuousTrajectory {
  using TimeTrajectory_t = std::array<Scalar, PredictLength + 1>;
  using StateTrajectory_t = std::array<Vector<Scalar, XDim>, PredictLength + 1>;
  using InputTrajectory_t = std::array<Vector<Scalar, UDim>, PredictLength>;
  /** 时间轨迹，尺寸 N+1。 */
  TimeTrajectory_t timeTrajectory;
  /** 状态向量轨迹，尺寸 N+1。 */
  StateTrajectory_t stateTrajectory;
  /** 输入向量轨迹，尺寸 N。 */
  InputTrajectory_t inputTrajectory;
};

/** 相加两条轨迹的状态和输入，不相加时间。 */
template <typename Scalar, int XDim, int UDim, size_t PredictLength>
ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> operator+(
    const ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>& lhs,
    const ContinuousTrajectory<Scalar, XDim, UDim, PredictLength>& rhs) {
  // 将 lhs 复制到 sum。
  ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> sum(lhs);

  for (size_t k = 0; k < sum.inputTrajectory.size(); ++k) {
    sum.inputTrajectory[k] += rhs.inputTrajectory[k];
  }

  // 状态求和。
  for (size_t k = 0; k < sum.stateTrajectory.size(); ++k) {
    sum.stateTrajectory[k] += rhs.stateTrajectory[k];
  }
  return sum;
}

/** 轨迹上一点的引用，不拥有状态-输入数据。
 */
template <typename Scalar, int XDim, int UDim>
struct TrajectoryRef {
  using StateVector_t = Vector<Scalar, XDim>;
  using InputVector_t = Vector<Scalar, UDim>;
  /** 时间 */
  Scalar t;
  /** 状态 */
  const StateVector_t& x;
  /** 输入。 */
  const InputVector_t& u;
};

/** 轨迹上一点状态的引用，不拥有状态
 * 数据。 */
template <typename Scalar, int XDim>
struct StateTrajectoryRef {
  using StateVector_t = Vector<Scalar, XDim>;
  /** 时间 */
  Scalar t;
  /** 状态 */
  const StateVector_t& x;
};

}  // namespace qp_solver
