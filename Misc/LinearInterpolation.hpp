#pragma once
#include "Types.hpp"
#include <array>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <cstddef>

namespace LinearInterpolation
{
  /**
   *  Find index into a sorted time Array
   *
   *  Indices are counted as follows:
   *          ------ | ----- | ---  ... ---    | -----
   *                t0     t1              t(n-1)
   *  Index     0        1      2   ...  (n-1)    n
   *
   *  Corner cases:
   *     - If time equal to a time in the timeArray is requested, the lower index is taken (e.g. t = t1 -> index = 1)
   *     - If multiple times in the timeArray are equal, the index before the first occurrence is taken.
   *       for example: if t1 = t2 = t3  and the requested time t <= t3 -> index = 1
   *
   *
   * @tparam Scalar : numerical type of time
   * @param timeArray : sorted time array to perform the lookup in
   * @param time : enquiry time
   * @return index between [0, size(timeArray)]
   */
  template <typename Scalar, size_t ArrayLength>
  size_t findIndexInTimeArray(const std::array<Scalar, ArrayLength> &timeArray, Scalar time)
  {
    size_t index = 0;
    for (index = 0; index < ArrayLength; ++index)
    {
      if (time < timeArray[index])
      {
        break;
      }
    }
    return index;
  }

  /**
   *  Find interval into a sorted time Array
   *
   *  Intervals are counted as follows:
   *           ------ | ----- | ---  ... ---    | -----
   *                 t0     t1              t(n-1)
   *  Interval  -1       0      1   ...  (n-2)    (n-1)
   *
   *  Corner cases are handled as in findIndexInTimeArray
   *
   * @tparam SCALAR : numerical type of time
   * @param timeArray : sorted time array to perform the lookup in
   * @param time : enquiry time
   * @return interval between [-1, size(timeArray)-1]
   */
  template <typename Scalar, size_t ArrayLength>
  int findIntervalInTimeArray(const std::array<Scalar, ArrayLength> &timeArray, Scalar time)
  {
    return findIndexInTimeArray(timeArray, time) - 1;
  }

  template <typename Scalar, size_t ArrayLength>
  std::pair<int, Scalar> timeSegment(Scalar enquiryTime, const std::array<Scalar, ArrayLength> &timeArray)
  {
    // corner cases (no time set OR single time element)
    if constexpr (ArrayLength <= 1)
    {
      return {0, Scalar(1.0)};
    }

    const int index = findIntervalInTimeArray(timeArray, enquiryTime);
    const int lastInterval = ArrayLength - 1;
    if (index >= 0)
    {
      if (index < lastInterval)
      {
        // interpolation : 0 <= index < lastInterval
        assert(enquiryTime <= timeArray[index + 1]); // assert upper bound of lookup
        assert(timeArray[index] <= enquiryTime);     // assert lower bound of lookup
        const Scalar intervalLength = timeArray[index + 1] - timeArray[index];
        const Scalar timeTillNext = timeArray[index + 1] - enquiryTime;

        // Normal case: interval is large enough for normal interpolation
        constexpr Scalar minIntervalTime = 2.0 * 1e-5;
        if (intervalLength > minIntervalTime)
        {
          return {index, (timeTillNext / intervalLength)};
        }

        // Take closes point for small time intervals
        if (timeTillNext < 0.5 * intervalLength)
        { // short interval, closest to time[index + 1]
          return {index, Scalar(0.0)};
        }
        else
        { // short interval, closest to time[index]
          return {index, Scalar(1.0)};
        }
      }
      else
      {
        // upper bound : index >= lastInterval
        return {std::max(lastInterval - 1, 0), Scalar(0.0)};
      }
    }
    else
    {
      // lower bound : index < 0
      return {0, Scalar(1.0)};
    }
  }

  template <typename Data, size_t ArrayLength>
  const Data &stdAccessFun(const std::array<Data, ArrayLength> &arr, size_t ind)
  {
    return arr[ind];
  }

  /******************************************************************************************************/
  /******************************************************************************************************/
  /******************************************************************************************************/
  template <typename Scalar, typename Data, size_t ArrayLength, class AccessFun>
  auto interpolate(const std::pair<int, Scalar> &indexAlpha, const std::array<Data, ArrayLength> &dataArray, AccessFun accessFun)
      -> std::decay_t<typename std::result_of<AccessFun(const std::array<Data, ArrayLength> &, size_t)>::type>
  {
    static_assert(ArrayLength > 0);

    if constexpr (ArrayLength > 1)
    {
      // Normal interpolation case
      int index = indexAlpha.first;
      Scalar alpha = indexAlpha.second;
      const auto &lhs = accessFun(dataArray, index);
      const auto &rhs = accessFun(dataArray, index + 1);

      return alpha * lhs + (Scalar(1.0) - alpha) * rhs;
    }
    else
    { // dataArray.size() == 1
      // Time vector has only 1 element -> Constant function
      return accessFun(dataArray, 0);
    }
  }

  template <typename Scalar, typename Data, size_t ArrayLength, class AccessFun>
  auto interpolate(const Scalar enquiryTime, const std::array<Scalar, ArrayLength> &timeArray, const std::array<Data, ArrayLength> &dataArray,
                   AccessFun accessFun) -> std::decay_t<typename std::result_of<AccessFun(const std::array<Data, ArrayLength> &, size_t)>::type>
  {
    return interpolate(timeSegment(enquiryTime, timeArray), dataArray, accessFun);
  }

  /******************************************************************************************************/
  /******************************************************************************************************/
  /******************************************************************************************************/
  template <typename Scalar, typename Data, size_t ArrayLength>
  Data interpolate(const std::pair<int, Scalar> &indexAlpha, const std::array<Data, ArrayLength> &dataArray)
  {
    return interpolate(indexAlpha, dataArray, stdAccessFun<Data, ArrayLength>);
  }

  template <typename Scalar, typename Data, size_t ArrayLength>
  Data interpolate(const Scalar enquiryTime, const std::array<Data, ArrayLength> &timeArray, const std::array<Data, ArrayLength> &dataArray)
  {
    return interpolate(enquiryTime, timeArray, dataArray, stdAccessFun<Data, ArrayLength>);
  }

}