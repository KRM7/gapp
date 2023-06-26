/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_BOUNDED_VALUE_HPP
#define GA_UTILITY_BOUNDED_VALUE_HPP

#include "utility.hpp"
#include <utility>
#include <limits>

namespace gapp::detail
{
    /* Class representing an arbitrary interval. */
    template<typename T>
    struct Interval
    {
        constexpr Interval(const T& left, const T& right, bool left_inclusive, bool right_inclusive) noexcept :
            left_(left), right_(right), left_inclusive_(left_inclusive), right_inclusive_(right_inclusive)
        {}

        constexpr bool contains(const T& value) const noexcept
        {
            if (value < left_  || (!left_inclusive_  && value == left_))  return false;
            if (value > right_ || (!right_inclusive_ && value == right_)) return false;
            return true;
        }

        T left_;
        T right_;
        bool left_inclusive_;
        bool right_inclusive_;
    };

    /* Class representing a value within an interval. */
    template<typename T, Interval<T> I>
    class BoundedValue
    {
    public:
        using value_type = T;

        constexpr /* implicit */ BoundedValue(value_type value) noexcept
        {
            GAPP_ASSERT(I.contains(value), "The value is outside of the allowed interval.");

            value_ = std::move(value);
        }

        constexpr /* implicit */ operator value_type() const noexcept { return value_; }

        static constexpr T lower_bound = I.left_;
        static constexpr T upper_bound = I.right_;
        static constexpr bool left_inclusive = I.left_inclusive_;
        static constexpr bool right_inclusive = I.right_inclusive_;

    private:
        value_type value_;
    };


} // namespace gapp::detail

namespace gapp
{
    /** Type representing values in the closed interval [Low, High]. */
    template<typename T, T Low, T High>
    using BoundedValue = detail::BoundedValue<T, detail::Interval<T>{ Low, High, true, true }>;

    /** Type representing values in the closed interval [0, Tmax]. */
    template<typename T>
    using NonNegative = detail::BoundedValue<T, detail::Interval<T>{ T(0), std::numeric_limits<T>::max(), true, true }>;

    /** Type representing values in the half-closed interval [Tmin, 0). */
    template<typename T>
    using Negative = detail::BoundedValue<T, detail::Interval<T>{ std::numeric_limits<T>::lowest(), T(0), true, false }>;

    /** Type representing values in the half-closed interval (0, Tmax]. */
    template<typename T>
    using Positive = detail::BoundedValue<T, detail::Interval<T>{ T(0), std::numeric_limits<T>::max(), false, true }>;


    /** Type representing a probability value in the closed interval [0.0, 1.0]. */
    using Probability = detail::BoundedValue<double, detail::Interval<double>{ 0.0, 1.0, true, true }>;

    constexpr Probability operator ""_p(long double arg) noexcept
    {
        return { double(arg) };
    }

    /** Type representing a value in the closed interval of [0.0, 1.0]. */
    template<typename T>
    using Normalized = detail::BoundedValue<T, detail::Interval<T>{ 0.0, 1.0, true, true }>;

} // namespace gapp

#endif // !GA_UTILITY_BOUNDED_VALUE_HPP