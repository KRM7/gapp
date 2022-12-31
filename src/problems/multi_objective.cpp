/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "multi_objective.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <numbers>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::problems
{
    using std::numbers::pi;

    {


    }


    std::vector<double> Kursawe::invoke(const std::vector<double>& vars) const
    {
        assert(!vars.empty());

        double f1 = 0.0;
        double f2 = 0.0;

        for (size_t i = 0; i < vars.size() - 1; i++)
        {
            f1 += -10.0 * std::exp(-0.2 * std::sqrt(std::pow(vars[i], 2) + std::pow(vars[i + 1], 2)));
            f2 += std::pow(std::abs(vars[i]), 0.8) + 5.0 * std::sin(std::pow(vars[i], 3));
        }
        f2 += std::pow(std::abs(vars.back()), 0.8) + 5.0 * std::sin(std::pow(vars.back(), 3));

        return { -f1, -f2 };
    }

    std::vector<double> ZDT1::invoke(const std::vector<double>& vars) const
    {
        assert(!vars.empty());

        double f1 = vars[0];

        double g = std::reduce(vars.begin() + 1, vars.end(), 0.0);
        g = 1.0 + 9.0 * g / (vars.size() - 1.0);

        double f2 = g - g * std::sqrt(f1 / g);

        return { -f1, -f2 };
    }

    std::vector<double> ZDT2::invoke(const std::vector<double>& vars) const
    {
        assert(!vars.empty());

        double f1 = vars[0];

        double g = std::reduce(vars.begin() + 1, vars.end(), 0.0);
        g = 1.0 + 9.0 * g / (vars.size() - 1.0);

        double f2 = g - f1 * f1 / g;

        return { -f1, -f2 };
    }

    std::vector<double> ZDT3::invoke(const std::vector<double>& vars) const
    {
        assert(!vars.empty());

        double f1 = vars[0];

        double g = std::reduce(vars.begin() + 1, vars.end(), 0.0);
        g = 1.0 + 9.0 * g / (vars.size() - 1.0);

        double f2 = g - g * std::sqrt(f1 / g) - f1 * std::sin(10.0 * pi * f1);

        return { -f1, -f2 };
    }

    std::vector<double> ZDT4::invoke(const std::vector<double>& vars) const
    {
        assert(!vars.empty());

        double f1 = vars[0];

        constexpr auto transform = [](double val) noexcept { return std::pow(val, 2) - 10.0 * std::cos(4.0 * pi * val); };

        double g = std::transform_reduce(vars.begin() + 1, vars.end(), 0.0, std::plus{}, transform);
        g = 1.0 + 10.0 * (vars.size() - 1.0) + g;

        double f2 = g - g * std::sqrt(f1 / g);

        return { -f1, -f2 };
    }

    std::vector<double> ZDT5::invoke(const std::vector<char>& vars) const
    {
        assert(vars.size() >= FIRST_BITS + REST_BITS);
        assert(vars.size() % REST_BITS == 0);

        double f1 = 1.0 + std::count(vars.begin(), vars.begin() + FIRST_BITS, BinaryGene{ 1 });

        double g = 0.0;
        for (auto first = vars.begin() + FIRST_BITS, last = first + REST_BITS; last < vars.end(); first += REST_BITS, last += REST_BITS)
        {
            size_t ones = std::count(first, last, BinaryGene{ 1 });
            g += (ones == REST_BITS) ? 1.0 : 2.0 + ones;
        }

        double f2 = g / f1;

        return { -f1, -f2 };
    }

    std::vector<double> ZDT6::invoke(const std::vector<double>& vars) const
    {
        assert(!vars.empty());

        double f1 = 1.0 - std::exp(-4.0 * vars[0]) * std::pow(std::sin(6.0 * pi * vars[0]), 6);

        double g = std::reduce(vars.begin() + 1, vars.end(), 0.0);
        g = 1.0 + 9.0 * std::pow(g / (vars.size() - 1.0), 0.25);

        double f2 = g - f1 * f1 / g;

        return { -f1, -f2 };
    }

} // namespace genetic_algorithm::problems