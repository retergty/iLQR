#pragma once
#include "Types.hpp"
#include <array>
#include "IntrusiveList.hpp"
#include <type_traits>
#include "LinearInterpolation.hpp"

template<typename Scalar>
struct Multiplier
{
  Multiplier() : Multiplier(0, 0) {}
  Multiplier(const Scalar penaltyArg, const Scalar lagrangianArg) : penalty(penaltyArg), lagrangian(lagrangianArg) {}

  Scalar penalty;
  Scalar lagrangian;
};

/**
 * The collection of Multiplier structure for all possible constraint terms in a particular time point.
 * stateEq : An array of state equality constraint terms.
 * stateIneq : An array of state inequality constraint terms.
 * stateInputEq : An array of state-input equality constraint terms.
 * stateInputIneq : An array of state-input inequality constraint terms.
 */
template<typename Scalar, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains>
struct MultiplierCollection
{
  // state equality
  std::array<Multiplier<Scalar>, StateEqConstrains> stateEq;
  // state inequality
  std::array<Multiplier<Scalar>, StateIneqConstrains> stateIneq;
  // state-input equality
  std::array<Multiplier<Scalar>, StateInputEqConstrains> stateInputEq;
  // state-input inequality
  std::array<Multiplier<Scalar>, StateInputIneqConstrains> stateInputIneq;

  /** Exchanges the values of MultiplierCollection. */
  void swap(MultiplierCollection& other) {
    stateEq.swap(other.stateEq);
    stateIneq.swap(other.stateIneq);
    stateInputEq.swap(other.stateInputEq);
    stateInputIneq.swap(other.stateInputIneq);
  }
};

namespace LinearInterpolation {

  /**
   * Linearly interpolates a trajectory of Multipliers.
   *
   * @param [in] indexAlpha : index and interpolation coefficient (alpha) pair.
   * @param [in] dataArray : A trajectory of MultiplierConstRef.
   * @return The interpolated Multiplier at indexAlpha.
   */
  template<typename Scalar, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains, size_t ArrayLen>
  Multiplier<Scalar> interpolate(const std::pair<int, Scalar>& indexAlpha, const std::array<Multiplier<Scalar>, ArrayLen>& dataArray)
  {
    const Scalar penalty = interpolate(
      indexAlpha, dataArray, [](const std::array<Multiplier<Scalar>, ArrayLen>& array, size_t t) -> const Scalar& { return array[t].penalty; });

    const Scalar lagrangian = interpolate(
      indexAlpha, dataArray, [](const std::array<Multiplier<Scalar>, ArrayLen> &array, size_t t) -> const Scalar& { return array[t].lagrangian; });

    return { penalty, lagrangian };
  }

  /**
   * Linearly interpolates a trajectory of MultiplierCollections.
   *
   * @param [in] indexAlpha : index and interpolation coefficient (alpha) pair.
   * @param [in] dataArray : A trajectory of MultiplierCollections.
   * @return The interpolated MultiplierCollection at indexAlpha.
   */
  template<typename Scalar, int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains, size_t ArrayLen>
  MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> interpolate(const std::pair<int, Scalar>& indexAlpha,
    const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& dataArray)
  {
    // number of terms
    MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains> out;

    // state equality
    for (size_t i = 0; i < StateEqConstrains; i++) {
      Scalar penalty = interpolate(indexAlpha, dataArray, [i](const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& array, size_t t) -> const Scalar& {
        return array[t].stateEq[i].penalty;
        });
      Scalar lagrangian = interpolate(indexAlpha, dataArray, [i](const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& array, size_t t) -> const Scalar& {
        return array[t].stateEq[i].lagrangian;
        });
      out.stateEq[i] = { penalty,lagrangian };
    }  // end of i loop

    // state inequality
    for (size_t i = 0; i < StateIneqConstrains; i++) {
      Scalar penalty = interpolate(indexAlpha, dataArray, [i](const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& array, size_t t) -> const Scalar& {
        return array[t].stateIneq[i].penalty;
        });
      Scalar lagrangian = interpolate(indexAlpha, dataArray, [i](const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& array, size_t t) -> const Scalar& {
        return array[t].stateIneq[i].lagrangian;
        });
      out.stateIneq[i] = { penalty,lagrangian };
    }  // end of i loop

    // state-input equality
    for (size_t i = 0; i < StateInputEqConstrains; i++) {
      Scalar penalty = interpolate(indexAlpha, dataArray, [i](const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& array, size_t t) -> const Scalar& {
        return array[t].stateInputEq[i].penalty;
        });
      Scalar lagrangian = interpolate(indexAlpha, dataArray, [i](const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& array, size_t t) -> const Scalar& {
        return array[t].stateInputEq[i].lagrangian;
        });
      out.stateInputEq[i] = { penalty,lagrangian };
    }  // end of i loop

    // state-input inequality
    for (size_t i = 0; i < StateInputIneqConstrains; i++) {
      Scalar penalty = interpolate(indexAlpha, dataArray, [i](const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& array, size_t t) -> const Scalar {
        return array[t].stateInputIneq[i].penalty;
        });
      Scalar lagrangian = interpolate(indexAlpha, dataArray, [i](const std::array<MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>, ArrayLen>& array, size_t t) -> const Scalar& {
        return array[t].stateInputIneq[i].lagrangian;
        });
      out.stateInputIneq[i] = { penalty,lagrangian };
    }  // end of i loop

    return out;
  }

}  // namespace LinearInterpolation