/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "many_objective.hpp"
#include "../utility/math.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <numbers>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::problems
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

    static constexpr double dtlz7_g(const_iterator first, const_iterator last) noexcept
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

        const auto theta = [g](double x) { return (g * x + 0.5) / (1.0 + g); };

        *fx.rbegin()   = std::sin(*first * pi / 2.0);
        *++fx.rbegin() = std::cos(*first * pi / 2.0);
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

        fx.back() = (1.0 + g) * fx.size();
        for (auto fxi = fx.begin(); fxi != fx.end() - 1; ++first, ++fxi)
        {
            fx.back() -= *first * (1.0 + std::sin(3.0 * pi * *first));
            *fxi = *first / (1.0 + g);
        }
        fx.back() /= (1.0 + g);

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


    DTLZ1::DTLZ1(size_t num_obj, size_t bits_per_var) :
        BenchmarkFunction("DTLZ1", num_obj + K - 1, num_obj, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_obj < 2) GA_THROW(std::invalid_argument, "The number of objectives must be at least 2.");

        ideal_point_ = math::Point(num_obj, 0.0);
        nadir_point_ = math::Point(num_obj, -0.5);

        optimum_ = std::vector(num_vars(), 0.5);
        std::fill(optimum_.begin(), optimum_.begin() + num_obj - 1, 0.0);
        optimal_value_ = math::Point(num_obj, 0.0);
        optimal_value_.back() = -0.5;
    }

    auto DTLZ1::invoke(const std::vector<RealGene>& vars) const -> FitnessVector
    {
        return dtlz1(vars, num_objectives());
    }


    DTLZ2::DTLZ2(size_t num_obj, size_t bits_per_var) :
        BenchmarkFunction("DTLZ2", num_obj + K - 1, num_obj, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_obj < 2) GA_THROW(std::invalid_argument, "The number of objectives must be at least 2.");

        ideal_point_ = math::Point(num_obj, 0.0);
        nadir_point_ = math::Point(num_obj, -1.0);

        optimum_ = std::vector(num_vars(), 0.5);
        std::fill(optimum_.begin(), optimum_.begin() + num_obj - 1, 0.0);
        optimal_value_ = math::Point(num_obj, 0.0);
        optimal_value_[0] = -1.0;
    }

    auto DTLZ2::invoke(const std::vector<RealGene>& vars) const -> FitnessVector
    {
        return dtlz2(vars, num_objectives());
    }


    DTLZ3::DTLZ3(size_t num_obj, size_t bits_per_var) :
        BenchmarkFunction("DTLZ3", num_obj + K - 1, num_obj, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_obj < 2) GA_THROW(std::invalid_argument, "The number of objectives must be at least 2.");

        ideal_point_ = math::Point(num_obj, 0.0);
        nadir_point_ = math::Point(num_obj, -1.0);

        optimum_ = std::vector(num_vars(), 0.5);
        std::fill(optimum_.begin(), optimum_.begin() + num_obj - 1, 0.0);
        optimal_value_ = math::Point(num_obj, 0.0);
        optimal_value_[0] = -1.0;
    }

    auto DTLZ3::invoke(const std::vector<RealGene>& vars) const -> FitnessVector
    {
        return dtlz3(vars, num_objectives());
    }


    DTLZ4::DTLZ4(size_t num_obj, size_t bits_per_var) :
        BenchmarkFunction("DTLZ4", num_obj + K - 1, num_obj, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_obj < 2) GA_THROW(std::invalid_argument, "The number of objectives must be at least 2.");

        ideal_point_ = math::Point(num_obj, 0.0);
        nadir_point_ = math::Point(num_obj, -1.0);

        optimum_ = std::vector(num_vars(), 0.5);
        optimal_value_ = math::Point(num_obj, 0.0);
        optimal_value_[0] = -1.0;
    }

    auto DTLZ4::invoke(const std::vector<RealGene>& vars) const -> FitnessVector
    {
        return dtlz4(vars, num_objectives());
    }


    DTLZ5::DTLZ5(size_t num_obj, size_t bits_per_var) :
        BenchmarkFunction("DTLZ5", num_obj + K - 1, num_obj, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_obj < 2) GA_THROW(std::invalid_argument, "The number of objectives must be at least 2.");

        ideal_point_ = math::Point(num_obj, 0.0);
        nadir_point_ = math::Point(num_obj, 0.0);
        for (size_t i = 0; i < num_obj; i++)
        {
            nadir_point_[i] = -1.0 / std::pow(std::sqrt(2), num_obj - 1 - i);
        }
        nadir_point_[0] = nadir_point_[1];

        optimum_ = std::vector(num_vars(), 0.5);
        std::fill(optimum_.begin(), optimum_.begin() + num_obj - 1, 0.0);
        optimal_value_ = nadir_point_;
        optimal_value_.back() = 0.0;
    }

    auto DTLZ5::invoke(const std::vector<RealGene>& vars) const -> FitnessVector
    {
        return dtlz5(vars, num_objectives());
    }


    DTLZ6::DTLZ6(size_t num_obj, size_t bits_per_var) :
        BenchmarkFunction("DTLZ6", num_obj + K - 1, num_obj, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_obj < 2) GA_THROW(std::invalid_argument, "The number of objectives must be at least 2.");

        ideal_point_ = math::Point(num_obj, 0.0);
        nadir_point_ = math::Point(num_obj, 0.0);
        for (size_t i = 0; i < num_obj; i++)
        {
            nadir_point_[i] = -1.0 / std::pow(std::sqrt(2), num_obj - 1 - i);
        }
        nadir_point_[0] = nadir_point_[1];

        optimum_ = std::vector(num_vars(), 0.0);
        optimal_value_ = nadir_point_;
        optimal_value_.back() = 0.0;
    }

    auto DTLZ6::invoke(const std::vector<RealGene>& vars) const -> FitnessVector
    {
        return dtlz6(vars, num_objectives());
    }


    DTLZ7::DTLZ7(size_t num_obj, size_t bits_per_var) :
        BenchmarkFunction("DTLZ7", num_obj + K - 1, num_obj, Bounds{ 0.0, 1.0 }, bits_per_var)
    {
        if (num_obj < 2) GA_THROW(std::invalid_argument, "The number of objectives must be at least 2.");

        optimum_ = std::vector(num_vars(), 0.0);
        optimal_value_ = math::Point(num_obj, 0.0);
        optimal_value_.back() = -2.0 * num_obj;

        ideal_point_ = math::Point(num_obj, 0.0);
        ideal_point_.back() = -0.307004 * num_obj - 1.692996;
        nadir_point_ = math::Point(num_obj, -1.0);
        nadir_point_.back() = -2.0 * num_obj;
    }

    auto DTLZ7::invoke(const std::vector<RealGene>& vars) const -> FitnessVector
    {
        return dtlz7(vars, num_objectives());
    }

} // namespace genetic_algorithm::problems