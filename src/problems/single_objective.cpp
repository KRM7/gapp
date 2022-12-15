/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "single_objective.hpp"
#include <vector>
#include <numeric>
#include <numbers>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::problems
{
    using std::numbers::pi;
    using std::numbers::e;

    std::vector<RealGene> BenchmarkFunctionReal1::convert(const std::vector<BinaryGene>& binary_chrom) const
    {
        assert((binary_chrom.size() / var_bits_) == num_vars());

        std::vector<RealGene> vars(num_vars());

        for (size_t i = 0; i < vars.size(); i++)
        {
            const auto first = binary_chrom.begin() + i * var_bits_;
            const auto last = binary_chrom.begin() + (i + 1) * var_bits_;

            const RealGene val = std::accumulate(first, last, 0.0, [](RealGene acc, BinaryGene bit) noexcept
            {
                return (acc * 2) + bit;
            });

            vars[i] = val / (std::pow(2.0, var_bits_) - 1); // use double to avoid integer overflow
            vars[i] *= bounds_[i].upper - bounds_[i].lower;
            vars[i] += bounds_[i].lower;
        }

        return vars;
    }


    std::vector<double> Sphere::invoke(const std::vector<double>& vars) const
    {
        return { -std::inner_product(vars.begin(), vars.end(), vars.begin(), 0.0) };
    }

    std::vector<double> Rastrigin::invoke(const std::vector<double>& vars) const
    {
        double fx = 10.0 * vars.size();
        for (double var : vars)
        {
            fx += std::pow(var, 2) - 10.0 * std::cos(2.0 * pi * var);
        }

        return { -fx };
    }

    std::vector<double> Rosenbrock::invoke(const std::vector<double>& vars) const
    {
        double fx = 0.0;
        for (size_t i = 0; i < vars.size() - 1; i++)
        {
            fx += 100.0 * std::pow((vars[i + 1] - std::pow(vars[i], 2)), 2) + std::pow(vars[i] - 1.0, 2);
        }

        return { -fx };
    }

    std::vector<double> Schwefel::invoke(const std::vector<double>& vars) const
    {
        double fx = 418.98288727215591 * vars.size();
        for (double var : vars)
        {
            fx -= var * std::sin(std::sqrt(std::abs(var)));
        }

        return { -fx };
    }

    std::vector<double> Griewank::invoke(const std::vector<double>& vars) const
    {
        double fx = 4000.0;
        double fm = 1.0;
        for (size_t i = 0; i < vars.size(); i++)
        {
            fx += std::pow(vars[i], 2);
            fm *= std::cos(vars[i] / std::sqrt(i + 1.0));
        }

        return { -fx / 4000.0 + fm };
    }

    std::vector<double> Ackley::invoke(const std::vector<double>& vars) const
    {
        double f1 = 0.0;
        double f2 = 0.0;
        for (double var : vars)
        {
            f1 += std::pow(var, 2);
            f2 += std::cos(2.0 * pi * var);
        }
        f1 = std::exp(-0.2 * std::sqrt(f1 / vars.size()));
        f2 = std::exp(f2 / vars.size());

        double fx = -20.0 * f1 - f2 + 20.0 + e;

        return { -fx };
    }

    std::vector<double> Levy::invoke(const std::vector<double>& vars) const
    {
        assert(!vars.empty());

        constexpr auto fw = [](double val) noexcept { return 0.75 + val / 4.0; };

        double fx = std::pow(std::sin(pi * fw(vars[0])), 2);
        for (size_t i = 0; i < vars.size() - 1; i++)
        {
            fx += std::pow(fw(vars[i]) - 1.0, 2) * (1.0 + 10.0 * std::pow(std::sin(pi * fw(vars[i]) + 1.0), 2));
        }
        fx += std::pow(fw(vars.back()) - 1.0, 2) * (1.0 + std::pow(std::sin(2.0 * pi * fw(vars.back())), 2));

        return { -fx };
    }

} // namespace genetic_algorithm::problems