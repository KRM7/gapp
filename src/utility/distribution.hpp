/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_DISTRIBUTION_HPP
#define GAPP_UTILITY_DISTRIBUTION_HPP

#include "bit.hpp"
#include <random>
#include <type_traits>
#include <iosfwd>
#include <cstdint>

namespace gapp::detail
{
    class uniform_bool_distribution
    {
    public:
        using result_type = bool;

        struct param_type
        {
            using distribution_type = uniform_bool_distribution;
            friend constexpr bool operator==(const param_type&, const param_type&) = default;
        };

        constexpr uniform_bool_distribution() noexcept = default;
        explicit constexpr uniform_bool_distribution(const param_type&) noexcept {}

        template<std::uniform_random_bit_generator Generator>
        [[nodiscard]] result_type operator()(Generator& generator) noexcept(std::is_nothrow_invocable_v<Generator>)
        {
            if (bit_pool_ == 1)
            {
                bit_pool_ = generator() | detail::msb_mask<typename Generator::result_type>;
            }

            const bool bit = detail::lsb(bit_pool_);
            bit_pool_ >>= 1;

            return bit;
        }

        template<std::uniform_random_bit_generator Generator>
        [[nodiscard]] result_type operator()(Generator& generator, const param_type&) noexcept(std::is_nothrow_invocable_v<Generator>)
        {
            return operator()(generator);
        }

        constexpr void reset() noexcept { bit_pool_ = 1; }

        static constexpr result_type min() noexcept { return false; }
        static constexpr result_type max() noexcept { return true; }

        static constexpr void param(const param_type&) noexcept {}
        static constexpr param_type param() noexcept { return {}; }

        friend constexpr bool operator==(const uniform_bool_distribution&, const uniform_bool_distribution&) = default;

        template<typename CharT, typename Traits>
        friend std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const uniform_bool_distribution& dist)
        {
            os << dist.bit_pool_;
        }

        template<typename CharT, typename Traits>
        friend std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits>& is, uniform_bool_distribution& dist)
        {
            is >> dist.bit_pool_;
        }

    private:
        std::uint64_t bit_pool_ = 1;
    };

} // namespace gapp::detail

#endif // !GAPP_UTILITY_DISTRIBUTION_HPP
