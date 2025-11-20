#pragma once

namespace numeric_traits
{
    /**
     * If a floating number v is in limit epsilon vicinity of w, then we assume that v approaches to w.
     * This limit_epsilon value should be greater than the weak_epsilon value.
     *
     * @return limit epsilon value.
     */
    template <typename T>
    constexpr T limitEpsilon()
    {
        return T(1e-6);
    }

    /**
     * Defines the precision during comparison.
     *
     * @return weak epsilon value
     */
    template <typename T>
    constexpr T weakEpsilon()
    {
        return T(1e-9);
    }

} // namespace numeric_traits