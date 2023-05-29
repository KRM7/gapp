/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "single_objective.hpp"
#include "../utility/utility.hpp"
#include "../utility/functional.hpp"
#include <vector>
#include <numeric>
#include <numbers>
#include <cmath>
#include <cstddef>

namespace genetic_algorithm::problems
{
    using std::numbers::pi;
    using std::numbers::e;

    auto Sphere::invoke(const std::vector<double>& vars) const -> FitnessVector
    {
        return { -std::inner_product(vars.begin(), vars.end(), vars.begin(), 0.0) };
    }

    auto Rastrigin::invoke(const std::vector<double>& vars) const -> FitnessVector
    {
        double fx = 10.0 * vars.size();
        for (double var : vars)
        {
            fx += std::pow(var, 2) - 10.0 * std::cos(2.0 * pi * var);
        }

        return { -fx };
    }

    auto Rosenbrock::invoke(const std::vector<double>& vars) const -> FitnessVector
    {
        GA_ASSERT(!vars.empty());

        double fx = 0.0;
        for (size_t i = 0; i < vars.size() - 1; i++)
        {
            fx += 100.0 * std::pow((vars[i + 1] - std::pow(vars[i], 2)), 2) + std::pow(vars[i] - 1.0, 2);
        }

        return { -fx };
    }

    auto Schwefel::invoke(const std::vector<double>& vars) const -> FitnessVector
    {
        double fx = 418.98288727215591 * vars.size();
        for (double var : vars)
        {
            fx -= var * std::sin(std::sqrt(std::abs(var)));
        }

        return { -fx };
    }

    auto Griewank::invoke(const std::vector<double>& vars) const -> FitnessVector
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

    auto Ackley::invoke(const std::vector<double>& vars) const -> FitnessVector
    {
        GA_ASSERT(!vars.empty());

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

    auto Levy::invoke(const std::vector<double>& vars) const -> FitnessVector
    {
        GA_ASSERT(!vars.empty());

        constexpr auto fw = detail::multiply_add(0.25, 0.75);

        double fx = std::pow(std::sin(pi * fw(vars[0])), 2);
        for (size_t i = 0; i < vars.size() - 1; i++)
        {
            fx += std::pow(fw(vars[i]) - 1.0, 2) * (1.0 + 10.0 * std::pow(std::sin(pi * fw(vars[i]) + 1.0), 2));
        }
        fx += std::pow(fw(vars.back()) - 1.0, 2) * (1.0 + std::pow(std::sin(2.0 * pi * fw(vars.back())), 2));

        return { -fx };
    }

} // namespace genetic_algorithm::problems