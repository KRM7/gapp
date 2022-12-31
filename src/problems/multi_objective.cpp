/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "multi_objective.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::problems
{
    using std::numbers::pi;

    Kursawe::Kursawe(size_t num_vars, size_t bits_per_var) :
        BenchmarkFunctionRealN("Kursawe", 2, num_vars, Bounds{ -5.0, 5.0 }, bits_per_var)
    {
        if (num_vars < 2) GA_THROW(std::invalid_argument, "Number of variables must be at least 2.");

        optimum_ = std::vector(num_vars, 0.0);
        optimal_value_ = { 10.0 * double(num_vars - 1), 0.0 };

        ideal_point_ = { 10.0 * double(num_vars - 1), 3.85 * double(num_vars - 1) + 4.0 };
        nadir_point_ = { 7.25 * double(num_vars - 1), 0.0};
    }

    std::vector<double> Kursawe::invoke(const std::vector<double>& vars) const
    {
        assert(!vars.empty());

        double f1 = 0.0;
        double f2 = 0.0;

        for (size_t i = 0; i < vars.size() - 1; i++)
        {
            f1 -= 10.0 * std::exp(-0.2 * std::sqrt(std::pow(vars[i], 2) + std::pow(vars[i + 1], 2)));
            f2 += std::pow(std::abs(vars[i]), 0.8) + 5.0 * std::sin(std::pow(vars[i], 3));
        }
        f2 += std::pow(std::abs(vars.back()), 0.8) + 5.0 * std::sin(std::pow(vars.back(), 3));

        return { -f1, -f2 };
    }


    ZDT1::ZDT1(size_t num_vars, size_t bits_per_var) :
        BenchmarkFunctionRealN("ZDT1", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_vars < 2) GA_THROW(std::invalid_argument, "Number of variables must be at least 2.");

        optimum_ = std::vector(num_vars, 0.0);
        optimal_value_ = { 0.0, -1.0 };

        ideal_point_ = {  0.0,  0.0 };
        nadir_point_ = { -1.0, -1.0 };
    }

    std::vector<double> ZDT1::invoke(const std::vector<double>& vars) const
    {
        assert(vars.size() > 1);

        double f1 = vars[0];

        double g = std::reduce(vars.begin() + 1, vars.end(), 0.0);
        g = 1.0 + 9.0 * g / (vars.size() - 1.0);

        double f2 = g - g * std::sqrt(f1 / g);

        return { -f1, -f2 };
    }


    ZDT2::ZDT2(size_t num_vars, size_t bits_per_var) :
        BenchmarkFunctionRealN("ZDT2", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_vars < 2) GA_THROW(std::invalid_argument, "Number of variables must be at least 2.");

        optimum_ = std::vector(num_vars, 0.0);
        optimal_value_ = { 0.0, -1.0 };

        ideal_point_ = {  0.0,  0.0 };
        nadir_point_ = { -1.0, -1.0 };
    }

    std::vector<double> ZDT2::invoke(const std::vector<double>& vars) const
    {
        assert(vars.size() > 1);

        double f1 = vars[0];

        double g = std::reduce(vars.begin() + 1, vars.end(), 0.0);
        g = 1.0 + 9.0 * g / (vars.size() - 1.0);

        double f2 = g - f1 * f1 / g;

        return { -f1, -f2 };
    }


    ZDT3::ZDT3(size_t num_vars, size_t bits_per_var) :
        BenchmarkFunctionRealN("ZDT3", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_vars < 2) GA_THROW(std::invalid_argument, "Number of variables must be at least 2.");

        optimum_ = std::vector(num_vars, 0.0);
        optimal_value_ = { 0.0, -1.0 };

        ideal_point_ = {  0.0,   0.8 };
        nadir_point_ = { -0.85, -1.0 };
    }

    std::vector<double> ZDT3::invoke(const std::vector<double>& vars) const
    {
        assert(vars.size() > 1);

        double f1 = vars[0];

        double g = std::reduce(vars.begin() + 1, vars.end(), 0.0);
        g = 1.0 + 9.0 * g / (vars.size() - 1.0);

        double f2 = g - g * std::sqrt(f1 / g) - f1 * std::sin(10.0 * pi * f1);

        return { -f1, -f2 };
    }


    ZDT4::ZDT4(size_t num_vars, size_t bits_per_var) :
        BenchmarkFunctionRealN("ZDT4", 2, num_vars, Bounds{ -5.0, 5.0 }, bits_per_var)
    {
        if (num_vars < 2) GA_THROW(std::invalid_argument, "Number of variables must be at least 2.");

        bounds_[0] = { 0.0, 1.0 };

        optimum_ = std::vector(num_vars, 0.0);
        optimal_value_ = { 0.0, -1.0 };

        ideal_point_ = {  0.0,  0.0 };
        nadir_point_ = { -1.0, -1.0 };
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


    ZDT5::ZDT5(size_t num_vars) :
        BenchmarkFunctionBinaryN("ZDT5", 2, FIRST_BITS + (num_vars - 1) * REST_BITS)
    {
        if (num_vars < 2) GA_THROW(std::invalid_argument, "Number of variables must be at least 2.");

        optimum_ = std::vector(this->num_vars(), Gene{ 1 });
        optimal_value_ = { -(FIRST_BITS + 1.0), -(num_vars - 1.0) / (FIRST_BITS + 1.0) };

        ideal_point_ = { -1.0, -(num_vars - 1.0) / (FIRST_BITS + 1.0) };
        nadir_point_ = { -(FIRST_BITS + 1.0), -(num_vars - 1.0) };
    }

    std::vector<double> ZDT5::invoke(const std::vector<char>& vars) const
    {
        assert(vars.size() >= FIRST_BITS);
        assert(vars.size() % REST_BITS == 0);

        double f1 = 1.0 + std::count(vars.begin(), vars.begin() + FIRST_BITS, Gene{ 1 });

        double g = 0.0;
        for (auto first = vars.begin() + FIRST_BITS; first <= vars.end() - REST_BITS; first += REST_BITS)
        {
            size_t ones = std::count(first, first + REST_BITS, Gene{ 1 });
            g += (ones == REST_BITS) ? 1.0 : 2.0 + ones;
        }

        double f2 = g / f1;

        return { -f1, -f2 };
    }


    ZDT6::ZDT6(size_t num_vars, size_t bits_per_var) :
        BenchmarkFunctionRealN("ZDT6", 2, num_vars, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_vars < 2) GA_THROW(std::invalid_argument, "Number of variables must be at least 2.");

        optimum_ = std::vector(num_vars, 0.0);
        optimal_value_ = { -1.0, 0.0 };

        ideal_point_ = {  0.0,  0.0 };
        nadir_point_ = { -1.0, -0.92 };
    }

    std::vector<double> ZDT6::invoke(const std::vector<double>& vars) const
    {
        assert(vars.size() > 1);

        double f1 = 1.0 - std::exp(-4.0 * vars[0]) * std::pow(std::sin(6.0 * pi * vars[0]), 6);

        double g = std::reduce(vars.begin() + 1, vars.end(), 0.0);
        g = 1.0 + 9.0 * std::pow(g / (vars.size() - 1.0), 0.25);

        double f2 = g - f1 * f1 / g;

        return { -f1, -f2 };
    }

} // namespace genetic_algorithm::problems