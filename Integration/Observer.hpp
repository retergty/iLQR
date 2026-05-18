/**
 * @file Observer.hpp
 * @brief 积分观测器：在积分过程中将 (时间, 状态) 写入给定缓冲区。
 */
#pragma once

#include "Types.hpp"

/**
 * @brief 观测器：将积分器输出的 (time, state)
 * 依次写入外部提供的时间与状态数组。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 */
template <typename Scalar, int XDim>
class Observer {
 public:
  /**
   * @brief 构造：绑定写入长度及状态/时间数组指针（可为空）。
   * @param [in] length 最大写入点数。
   * @param [in] stateTrajectoryPtr 状态轨迹缓冲区指针。
   * @param [in] timeTrajectoryPtr 时间轨迹缓冲区指针。
   */
  explicit Observer(int length,
                    Vector<Scalar, XDim>* stateTrajectoryPtr = nullptr,
                    Scalar* timeTrajectoryPtr = nullptr)
      : stateTrajectoryPtr_(stateTrajectoryPtr),
        timeTrajectoryPtr_(timeTrajectoryPtr),
        length_(length) {};

  /** @brief 析构函数。 */
  ~Observer() = default;

  /**
   * @brief 记录当前时间与状态到缓冲区（由积分器在每步调用）。
   * @param [in] state 当前状态。
   * @param [in] time 当前时间。
   */
  void observe(const Vector<Scalar, XDim>& state, const Scalar time) {
    if (now_ >= length_) return;
    if (timeTrajectoryPtr_ != nullptr) timeTrajectoryPtr_[now_] = time;
    if (stateTrajectoryPtr_ != nullptr) stateTrajectoryPtr_[now_] = state;
    now_++;
    return;
  }

  /** @brief 重置写入计数为 0。 */
  void clear() { now_ = 0; }

  /** @brief 返回已写入的点数。 */
  int getCount() const { return now_; }

 private:
  Vector<Scalar, XDim>* stateTrajectoryPtr_;
  Scalar* timeTrajectoryPtr_;
  int length_;
  int now_{0};
};
