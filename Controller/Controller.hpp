#pragma once
#include "Types.hpp"
#include <array>
#include "LinearInterpolation.hpp"

/**
 * Enum class for specifying controller type
 */
enum class ControllerType
{
  UNKNOWN,
  FEEDFORWARD,
  LINEAR,
  ONNX,
  BEHAVIORAL
};

/**
 * The base class for all controllers.
 */
template <typename Scalar, int XDimisions, int UDimisions>
class ControllerBase
{
public:
  /** Constructor */
  ControllerBase() = default;

  /** Default destructor. */
  virtual ~ControllerBase() = default;

  /**
   * @brief Computes the control command at a given time and state.
   *
   * @param [in] t: Current time.
   * @param [in] x: Current state.
   * @return Current input.
   */
  virtual Vector<Scalar, UDimisions> computeInput(Scalar t, const Vector<Scalar, XDimisions>& x) = 0;

  virtual Vector<Scalar, UDimisions> computeInput(size_t time_index, const Vector<Scalar, XDimisions>& x) = 0;
  /**
   * @brief Prints the type of controller
   * @return ControllerType: what type of controller this is
   */
  virtual ControllerType getType() const = 0;

  /**
   * @brief clears and reverts back to an empty controller.
   * Therefore, if empty() method is called, it will return true.
   */
  virtual void clear() = 0;

  /**
   * Returns whether the class contains any information.
   *
   * @return true if it contains no information, false otherwise.
   */
  virtual bool empty() const = 0;

  /**
   * Displays controller's data.
   */
  virtual void display() const {}
};
