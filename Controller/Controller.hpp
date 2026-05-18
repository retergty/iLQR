/**
 * @file Controller.hpp
 * @brief 控制器基类：按时间与状态计算控制输入，支持清空与类型查询。
 */
#pragma once
#include <array>

#include "LinearInterpolation.hpp"
#include "Types.hpp"

/** @brief 控制器类型枚举。 */
enum class ControllerType { UNKNOWN, FEEDFORWARD, LINEAR, ONNX, BEHAVIORAL };

/**
 * @brief 控制器基类：根据时间 t 与状态 x 计算控制 u，提供 clear/empty/display
 * 等接口。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class ControllerBase {
 public:
  /** @brief 默认构造。 */
  ControllerBase() = default;

  /** @brief 虚析构函数。 */
  virtual ~ControllerBase() = default;

  /**
   * @brief 在给定时间和状态下计算控制量。
   * @param [in] t 当前时间。
   * @param [in] x 当前状态。
   * @return 当前控制输入。
   */
  virtual Vector<Scalar, UDim> computeInput(
      Scalar t, const Vector<Scalar, XDim>& x) const = 0;

  /** @brief 按离散时间索引与状态计算控制（用于固定时间网格）。 */
  virtual Vector<Scalar, UDim> computeInput(
      size_t time_index, const Vector<Scalar, XDim>& x) const = 0;
  /**
   * @brief 返回控制器类型。
   * @return 控制器类型枚举值。
   */
  virtual ControllerType getType() const = 0;

  /**
   * @brief 清空并恢复为空控制器；之后 empty() 返回 true。
   */
  virtual void clear() = 0;

  /**
   * @brief 判断控制器是否不含任何有效信息。
   * @return 无信息返回 true，否则返回 false。
   */
  virtual bool empty() const = 0;

  /** @brief 可选：在控制台输出控制器数据，默认空实现。 */
  virtual void display() const {}
};
