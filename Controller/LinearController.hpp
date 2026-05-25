/**
 * @file LinearController.hpp
 * @brief 线性控制器：u(t,x) = K(t)*x + uff(t)，支持时间/状态插值。
 */
#pragma once
#include "Controller.hpp"
#include "LinearInterpolation.hpp"

/**
 * @brief 线性控制器，形式为 u[x,t] = k[t]*x + uff[t]，时间戳与增益/偏置为数组。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 * @tparam ArrayLen 时间节点数（时间戳与增益/偏置数组长度）。
 */
template <typename Scalar, int XDim, int UDim, size_t ArrayLen>
class LinearController final : public ControllerBase<Scalar, XDim, UDim> {
 public:
  /** @brief 默认构造，成员未初始化。 */
  LinearController() = default;

  /**
   * @brief 用给定时间戳、前馈偏置与反馈增益数组初始化控制器。
   * @param [in] controllerTime 时间戳数组。
   * @param [in] controllerBias 前馈偏置数组 uff。
   * @param [in] controllerGain 反馈增益数组 K。
   */
  LinearController(
      const std::array<Scalar, ArrayLen>& controllerTime,
      const std::array<Vector<Scalar, UDim>, ArrayLen>& controllerBias,
      const std::array<Matrix<Scalar, UDim, XDim>, ArrayLen>& controllerGain)
      : timeStamp_(controllerTime),
        biasArray_(controllerBias),
        gainArray_(controllerGain) {}

  /** @brief 拷贝构造。 */
  LinearController(const LinearController& other)
      : LinearController(other.timeStamp_, other.biasArray_, other.gainArray_) {
    deltaBiasArray_ = other.deltaBiasArray_;
  }

  /** @brief 拷贝赋值。 */
  LinearController& operator=(const LinearController& rhs) {
    timeStamp_ = rhs.timeStamp_;
    biasArray_ = rhs.biasArray_;
    gainArray_ = rhs.gainArray_;
    deltaBiasArray_ = rhs.deltaBiasArray_;

    return *this;
  }

  /** @brief 析构函数。 */
  ~LinearController() override = default;

  /**
   * @brief 设置控制律：时间戳、前馈偏置与反馈增益。
   * @param [in] controllerTime 时间戳数组。
   * @param [in] controllerBias 前馈偏置数组。
   * @param [in] controllerGain 反馈增益数组。
   */
  void setController(
      const std::array<Scalar, ArrayLen>& controllerTime,
      const std::array<Vector<Scalar, UDim>, ArrayLen>& controllerBias,
      const std::array<Matrix<Scalar, UDim, XDim>, ArrayLen>& controllerGain) {
    timeStamp_ = controllerTime;
    biasArray_ = controllerBias;
    gainArray_ = controllerGain;
  }

  /** @brief 按时间 t 与状态 x 计算控制：先时间插值得到 K 与 uff，再 u = uff +
   * K*x。 */
  Vector<Scalar, UDim> computeInput(
      Scalar t, const Vector<Scalar, XDim>& x) const override {
    const std::pair<int, Scalar> indexAlpha =
        LinearInterpolation::timeSegment(t, timeStamp_);

    Vector<Scalar, UDim> uff =
        LinearInterpolation::interpolate(indexAlpha, biasArray_);
    const Matrix<Scalar, UDim, XDim> k =
        LinearInterpolation::interpolate(indexAlpha, gainArray_);

    uff.noalias() += k * x;
    return uff;
  }

  /** @brief 返回控制器类型 LINEAR。 */
  ControllerType getType() const override { return ControllerType::LINEAR; }

  /** @brief 清空控制器：时间戳置零，增益与偏置置零；之后 empty() 返回 true。 */
  void clear() override {
    for (size_t i = 0; i < ArrayLen; ++i) {
      timeStamp_[i] = 0;
      biasArray_[i].setZero();
      deltaBiasArray_[i].setZero();
      gainArray_[i].setZero();
    }
  }

  /** @brief 判断是否为空：若所有时间戳近似为 0 则视为空。 */
  bool empty() const override {
    for (size_t i = 0; i < ArrayLen; ++i) {
      if (abs(timeStamp_[i]) > 1e-6f) {
        return false;
      }
    }
    return true;
  }

  /** @brief 返回时间节点数（数组长度）。 */
  constexpr static size_t size() { return ArrayLen; }

  /** @brief 与另一线性控制器交换时间戳、偏置、deltaBias、增益。 */
  void swap(LinearController& other) {
    timeStamp_.swap(other.timeStamp_);
    biasArray_.swap(other.biasArray_);
    deltaBiasArray_.swap(other.deltaBiasArray_);
    gainArray_.swap(other.gainArray_);
  }

 public:
  std::array<Scalar, ArrayLen> timeStamp_;
  std::array<Vector<Scalar, UDim>, ArrayLen> biasArray_;
  std::array<Vector<Scalar, UDim>, ArrayLen> deltaBiasArray_;
  std::array<Matrix<Scalar, UDim, XDim>, ArrayLen> gainArray_;
};
