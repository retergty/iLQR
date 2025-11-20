#pragma once
#include "Controller.hpp"
#include "LinearInterpolation.hpp"
/**
 * LinearController implements a time and state dependent controller of the
 * form u[x,t] = k[t] * x + uff[t]
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t ArrayLen>
class LinearController final : public ControllerBase<Scalar, XDimisions, UDimisions>
{
public:
  /** Constructor, leaves object uninitialized */
  LinearController() = default;

  /**
   * @brief Constructor initializes all required members of the controller.
   *
   * @param [in] controllerTime: Time stamp array of the controller
   * @param [in] controllerBias: The bias array.
   * @param [in] controllerGain: The feedback gain array.
   */
  LinearController(
      const std::array<Scalar, ArrayLen> &controllerTime,
      const std::array<Vector<Scalar, UDimisions>, ArrayLen> &controllerBias,
      const std::array<Matrix<Scalar, UDimisions, XDimisions>, ArrayLen> &controllerGain)
      : timeStamp_(controllerTime), biasArray_(controllerBias), gainArray_(controllerGain)
  {
  }

  /** Copy constructor */
  LinearController(const LinearController &other) : LinearController(other.timeStamp_, other.biasArray_, other.gainArray_)
  {
    deltaBiasArray_ = other.deltaBiasArray_;
  }

  /** Copy assignment (copy and swap idiom) */
  LinearController &operator=(const LinearController &rhs)
  {
    timeStamp_ = rhs.timeStamp_;
    biasArray_ = rhs.biasArray_;
    gainArray_ = rhs.gainArray_;
    deltaBiasArray_ = rhs.deltaBiasArray_;

    return *this;
  }

  /** Destructor */
  ~LinearController() override = default;

  /**
   * @brief setController Assign control law
   * @param [in] controllerTime: Time stamp array of the controller
   * @param [in] controllerBias: The bias array.
   * @param [in] controllerGain: The feedback gain array.
   */
  void setController(
      const std::array<Scalar, ArrayLen> &controllerTime,
      const std::array<Vector<Scalar, UDimisions>, ArrayLen> &controllerBias,
      const std::array<Matrix<Scalar, UDimisions, XDimisions>, ArrayLen> &controllerGain)
  {
    timeStamp_ = controllerTime;
    biasArray_ = controllerBias;
    gainArray_ = controllerGain;
  }

  Vector<Scalar, UDimisions> computeInput(Scalar t, const Vector<Scalar, XDimisions> &x) override
  {
    const std::pair<int, Scalar> indexAlpha = LinearInterpolation::timeSegment(t, timeStamp_);

    Vector<Scalar, UDimisions> uff = LinearInterpolation::interpolate(indexAlpha, biasArray_);
    const Matrix<Scalar, UDimisions, XDimisions> k = LinearInterpolation::interpolate(indexAlpha, gainArray_);

    uff.noalias() += k * x;
    return uff;
  }

  Vector<Scalar, UDimisions> computeInput(size_t time_index, const Vector<Scalar, XDimisions> &x) override
  {
    assert(time_index < ArrayLen);
    Vector<Scalar, UDimisions> uff = biasArray_[time_index];
    const Matrix<Scalar, UDimisions, XDimisions> &k = gainArray_[time_index];

    uff += k * x;
    return uff;
  }

  ControllerType getType() const override
  {
    return ControllerType::LINEAR;
  }

  /**
   * @brief clears and reverts back to an empty controller.
   * Therefore, if empty() method is called, it will return true.
   */
  void clear() override
  {
    for (size_t i = 0; i < ArrayLen; ++i)
    {
      timeStamp_[i] = 0;
      biasArray_[i].setZero();
      deltaBiasArray_[i].setZero();
      gainArray_[i].setZero();
    }
  }

  /**
   * Returns whether the class contains any information.
   *
   * @return true if it contains no information, false otherwise.
   */
  bool empty() const override
  {
    bool is_empty = true;
    for (size_t i = 0; i < ArrayLen; ++i)
    {
      if (abs(timeStamp_[i]) > 1e-6f)
      {
        is_empty = false;
      }
    }
    return is_empty;
  }

  constexpr size_t size() const
  {
    return ArrayLen;
  }
public:
  std::array<Scalar, ArrayLen> timeStamp_;
  std::array<Vector<Scalar, UDimisions>, ArrayLen> biasArray_;
  std::array<Vector<Scalar, UDimisions>, ArrayLen> deltaBiasArray_;
  std::array<Matrix<Scalar, UDimisions, XDimisions>, ArrayLen> gainArray_;
};
