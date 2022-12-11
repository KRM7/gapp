/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "many_objective.hpp"
#include <vector>
#include <numeric>
#include <functional>
#include <iterator>
#include <numbers>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::benchmark
{
    using std::numbers::pi;
    using const_iterator = std::vector<double>::const_iterator;


    /* DTLZ SUITE G FUNCTIONS */

    static inline double dtlz1_g(const_iterator first, const_iterator last) noexcept
    {
        assert(std::distance(first, last) > 0);

        double g = double(last - first);
        for (; first != last; ++first)
        {
            g += std::pow(*first - 0.5, 2) - std::cos(20.0 * pi * (*first - 0.5));
        }

        return 100.0 * g;
    }

    static inline double dtlz2_g(const_iterator first, const_iterator last) noexcept
    {
        assert(std::distance(first, last) > 0);

        double g = 0.0;
        for (; first != last; ++first)
        {
            g += std::pow(*first - 0.5, 2);
        }

        return g;
    }

    static constexpr auto dtlz3_g = dtlz1_g;

    static constexpr auto dtlz4_g = dtlz2_g;

    static constexpr auto dtlz5_g = dtlz2_g;

    static inline double dtlz6_g(const_iterator first, const_iterator last) noexcept
    {
        assert(std::distance(first, last) > 0);

        double g = 0.0;
        for (; first != last; ++first)
        {
            g += std::pow(*first, 0.1);
        }

        return g;
    }

    static inline double dtlz7_g(const_iterator first, const_iterator last) noexcept
    {
        assert(std::distance(first, last) > 0);

        return 1.0 + 9.0 / double(last - first) * std::reduce(first, last, 0.0);
    }


    /* DTLZ SUITE F FUNCTIONS */

    static inline std::vector<double> dtlz1_f(const_iterator first, const_iterator last, double)
    {
        assert(std::distance(first, last) > 0);

        std::vector fx(last - first + 1, 0.5);

        for (auto fxi = fx.rbegin(); fxi != fx.rend() - 1; ++first, ++fxi)
        {
            *std::next(fxi) = *fxi * *first;
            *fxi *= 1.0 - *first;
        }

        return fx;
    }

    static inline std::vector<double> dtlz2_f(const_iterator first, const_iterator last, double)
    {
        assert(std::distance(first, last) > 0);

        std::vector fx(last - first + 1, 1.0);

        for (auto fxi = fx.rbegin(); fxi != fx.rend() - 1; ++first, ++fxi)
        {
            *std::next(fxi) = *fxi * std::cos(*first * pi / 2.0);
            *fxi *= std::sin(*first * pi / 2.0);
        }

        return fx;
    }

    static constexpr auto dtlz3_f = dtlz2_f;

    static inline std::vector<double> dtlz4_f(const_iterator first, const_iterator last, double)
    {
        assert(std::distance(first, last) > 0);

        std::vector fx(last - first + 1, 1.0);

        for (auto fxi = fx.rbegin(); fxi != fx.rend() - 1; ++first, ++fxi)
        {
            *std::next(fxi) = *fxi * std::cos(std::pow(*first, 100) * pi / 2.0);
            *fxi *= std::sin(std::pow(*first, 100) * pi / 2.0);
        }

        return fx;
    }

    static inline std::vector<double> dtlz5_f(const_iterator first, const_iterator last, double g)
    {
        assert(std::distance(first, last) > 0);

        std::vector fx(last - first + 1, 1.0);

        const auto theta = [&](double x) { return (1.0 + 2.0 * g * x) / (2.0 * (1.0 + g)); };

        *std::next(fx.rbegin()) = *fx.rbegin() * std::cos(*first * pi / 2.0);
        *fx.rbegin() *= std::sin(*first * pi / 2.0);
        ++first;

        for (auto fxi = fx.rbegin() + 1; fxi != fx.rend() - 1; ++first, ++fxi)
        {
            *std::next(fxi) = *fxi * std::cos(theta(*first) * pi / 2.0);
            *fxi *= std::sin(theta(*first) * pi / 2.0);
        }

        return fx;
    }

    static constexpr auto dtlz6_f = dtlz5_f;

    static inline std::vector<double> dtlz7_f(const_iterator first, const_iterator last, double g)
    {
        assert(std::distance(first, last) > 0);

        std::vector fx(last - first + 1, 0.0);

        fx.back() = double(fx.size());
        for (auto fxi = fx.begin(); fxi != fx.end() - 1; ++first, ++fxi)
        {
            fx.back() -= *first / (1.0 + g) * (1.0 + std::sin(3.0 * pi * *first));
            *fxi = *first / (1.0 + g);
        }

        return fx;
    }


    /* DTLZ SUITE FUNCTIONS */

    template<auto F, auto G>
    static inline std::vector<double> dtlz(const std::vector<double>& vars, size_t num_obj)
    {
        const auto middle = vars.begin() + num_obj - 1;

        const double g = G(middle, vars.end());
        std::vector fx = F(vars.begin(), middle, g);

        /* Maximization. */
        for (double& val : fx) { val *= -(1.0 + g); }

        return fx;
    }

    static constexpr auto dtlz1 = dtlz<dtlz1_f, dtlz1_g>;
    static constexpr auto dtlz2 = dtlz<dtlz2_f, dtlz2_g>;
    static constexpr auto dtlz3 = dtlz<dtlz3_f, dtlz3_g>;
    static constexpr auto dtlz4 = dtlz<dtlz4_f, dtlz4_g>;
    static constexpr auto dtlz5 = dtlz<dtlz5_f, dtlz5_g>;
    static constexpr auto dtlz6 = dtlz<dtlz6_f, dtlz6_g>;
    static constexpr auto dtlz7 = dtlz<dtlz7_f, dtlz7_g>;


    std::vector<double> DTLZ1::invoke(const std::vector<double>& vars) const
    {
        return dtlz1(vars, num_obj());
    }

    std::vector<double> DTLZ2::invoke(const std::vector<double>& vars) const
    {
        return dtlz2(vars, num_obj());
    }

    std::vector<double> DTLZ3::invoke(const std::vector<double>& vars) const
    {
        return dtlz3(vars, num_obj());
    }

    std::vector<double> DTLZ4::invoke(const std::vector<double>& vars) const
    {
        return dtlz4(vars, num_obj());
    }

    std::vector<double> DTLZ5::invoke(const std::vector<double>& vars) const
    {
        return dtlz5(vars, num_obj());
    }

    std::vector<double> DTLZ6::invoke(const std::vector<double>& vars) const
    {
        return dtlz6(vars, num_obj());
    }

    std::vector<double> DTLZ7::invoke(const std::vector<double>& vars) const
    {
        return dtlz7(vars, num_obj());
    }

} // namespace genetic_algorithm::benchmark