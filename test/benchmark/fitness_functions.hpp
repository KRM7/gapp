/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

/*
* Fitness functions for testing the genetic algorithms.
* Includes some functions for both single- and multi-objective algorithms.
* All objective functions are to be maximized using the GAs.
* Note: All of the functions are modified for maximization where needed, and return vector<double>.
*/

#ifndef FITNESS_FUNCTIONS_HPP
#define FITNESS_FUNCTIONS_HPP

#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <numeric>
#include <algorithm>
#include <cassert>
#include "../src/genetic_algorithm.hpp"
#include "benchmark_utils.hpp"

using GeneBounds = genetic_algorithm::GA<double>::GeneBounds;

constexpr double PI = 3.14159265358979323846; // TODO include numbers

/* Single-objective fitness functions. */

/*
* Implementation of the Rastrigin function for any number of dimensions.
* Evaluated on x_i = [-5.12, 5.12].
* The global optimum of the function is f(x) = 0, at x = (0, 0, ... , 0).
*/
class Rastrigin
{
public:

    explicit Rastrigin(size_t nvars = 10) noexcept : // TODO base ctors that take all params
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const // TODO asserts in base
    {
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double fx = 10.0 * x.size();
        for (size_t i = 0; i < x.size(); i++)
        {
            fx += std::pow(x[i], 2) - 10.0 * cos(2 * PI * x[i]);
        }

        return { -fx };	/* For maximization. */
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const // TODO base
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 10; // TODO base
    size_t var_bits = 32; // TODO base + nobj

    constexpr static size_t num_obj() noexcept { return 1; } // TODO base with var
    constexpr static double lbound() noexcept { return -5.12; }
    constexpr static double ubound() noexcept { return 5.12; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); } // TODO in base class
    constexpr static double optimal_value() noexcept { return 0.0; } // TODO return vector
    constexpr static double optimal_x() noexcept { return 0.0; } // TODO return vector
    constexpr static double intval() noexcept { return ubound() - lbound(); }  // TODO in base class
};

/*
* Implementation of the Rosenbrock function for any dimensions.
* Evaluated on x_i = [-2.048, 2.048].
* The global optimum of the function is f(x) = 0, at x = (1, 1, ... , 1).
*/
class Rosenbrock
{
public:

    explicit Rosenbrock(size_t nvars = 3) noexcept :
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double fx = 0.0;
        for (size_t i = 0; i < x.size() - 1; i++)
        {
            fx += 100.0 * pow((x[i + 1] - pow(x[i], 2)), 2) + pow(1 - x[i], 2);
        }

        return { -fx }; /* For maximization. */
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars =  3;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -2.048; }
    constexpr static double ubound() noexcept { return 2.048; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double optimal_value() noexcept { return 0.0; }
    constexpr static double optimal_x() noexcept { return 1.0; }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};

/*
* Implementation of the Schwefel function for any number of dimensions.
* Evaluated on x_i = [-500, 500].
* The global optimum of the function is f(x) = 0, at x = (420.9687, 420.9687, ... , 420.9687).
*/
class Schwefel
{
public:

    explicit Schwefel(size_t nvars = 10) noexcept :
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double fx = 418.9829 * x.size();
        for (size_t i = 0; i < x.size(); i++)
        {
            fx -= x[i] * std::sin(std::sqrt(std::abs(x[i])));
        }

        return { -fx }; /* For maximization. */
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -500.0; }
    constexpr static double ubound() noexcept { return 500.0; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double optimal_value() noexcept { return 0.0; }
    constexpr static double optimal_x() noexcept { return 420.9687; }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};

/*
* Implementation of the Griewank function for any number of dimensions.
* Evaluated on x_i = [-600, 600].
* The global optimum of the function is f(x) = 0, at x = (0, 0, ... , 0).
*/
class Griewank
{
public:

    explicit Griewank(size_t nvars = 10) noexcept :
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double fx = 1.0;
        double fminus = 1.0;
        for (size_t i = 0; i < x.size(); i++)
        {
            fx += x[i] * x[i] / 4000.0;
            fminus *= std::cos(x[i] / std::sqrt(i + 1.0));
        }
        fx -= fminus;

        return { -fx }; /* For maximization. */
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -600.0; }
    constexpr static double ubound() noexcept { return 600.0; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double optimal_value() noexcept { return 0.0; }
    constexpr static double optimal_x() noexcept { return 0.0; }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};

/*
* Implementation of the Ackley function for any number of dimensions.
* Evaluated on x_i = [-32.768, 32.768].
* The global optimum of the function is f(x) = 0, at x = (0, 0, ... , 0).
*/
class Ackley
{
public:

    explicit Ackley(size_t nvars = 10) noexcept :
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = 0.0;
        double f2 = 0.0;
        for (size_t i = 0; i < x.size(); i++)
        {
            f1 += std::pow(x[i], 2);
            f2 += std::cos(2.0 * PI * x[i]);
        }
        f1 = std::exp(-0.2 * std::sqrt(f1 / num_vars));
        f2 = std::exp(f2 / num_vars);
        
        double fx = -20.0 * f1 - f2 + 20.0 + std::exp(1.0);

        return { -fx }; /* For maximization. */
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -32.768; }
    constexpr static double ubound() noexcept { return 32.768; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double optimal_value() noexcept { return 0.0; }
    constexpr static double optimal_x() noexcept { return 0.0; }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};


/* Multi-objective fitness functions. */

/*
* Implementation of the Kursawe function for any number of dimensions.
* Evaluated on x_i = [-5.0, 5.0].
* Multiobjective fitness function (objectives = 2).
*/
class KUR
{
public:

    explicit KUR(size_t nvars = 3) noexcept :
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() > 1);
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = 0.0;
        for (size_t i = 0; i < x.size() - 1; i++)
        {
            f1 += -10.0 * std::exp(-0.2 * std::sqrt(std::pow(x[i], 2) + std::pow(x[i + 1], 2)));
        }
        double f2 = 0.0;
        for (size_t i = 0; i < x.size(); i++)
        {
            f2 += std::pow(std::abs(x[i]), 0.8) + 5 * std::sin(std::pow(x[i], 3));
        }

        return { -f1, -f2 };
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 3;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 2; }
    constexpr static double lbound() noexcept { return -5.0; }
    constexpr static double ubound() noexcept { return 5.0; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};

/*
* Implementation of the ZDT2 function for any number of dimensions.
* Evaluated on x_i = [0.0, 1.0].
* Multiobjective fitness function (objectives = 2).
*/
class ZDT2
{
public:

    explicit ZDT2(size_t nvars = 30) noexcept :
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() > 1);
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = x[0];
        double g = std::accumulate(x.begin() + 1, x.end(), 0.0);
        g = 9 * g / (x.size() - 1.0) + 1;
        double f2 = g * (1.0 - pow(x[0] / g, 2));

        return { -f1, -f2 };
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 30;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 2; }
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};

/*
* Implementation of the ZDT3 function for any number of dimensions.
* Evaluated on x_i = [0.0, 1.0].
* Multiobjective fitness function (objectives = 2).
*/
class ZDT3
{
public:

    explicit ZDT3(size_t nvars = 30) noexcept :
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() > 1);
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = x[0];
        double g = std::accumulate(x.begin() + 1, x.end(), 0.0);
        g = 1 + 9 * g / (x.size() - 1.0);
        double f2 = g * (1 - std::sqrt(x[0] / g) - (x[0] / g) * std::sin(10 * PI * x[0]));

        return { -f1, -f2 };
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 30;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 2; }
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};

/*
* Implementation of the ZDT6 function for any number of dimensions.
* Evaluated on x_i = [0.0, 1.0].
* Multiobjective fitness function (objectives = 2).
*/
class ZDT6
{
public:

    explicit ZDT6(size_t nvars = 10) noexcept :
        num_vars(nvars) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() > 1);
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = 1 - std::exp(-4 * x[0]) * std::pow(std::sin(6 * PI * x[0]), 6);
        double g = std::accumulate(x.begin() + 1, x.end(), 0.0);
        g = 9.0 * std::pow(g / (x.size() - 1.0), 0.25) + 1;
        double f2 = g * (1 - std::pow(f1 / g, 2));

        return { -f1, -f2 };
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 2; }
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};


/* Many-objective fitness functions. */

/* 
* Implementation of the DTLZ1 function for any number of dimensions and objectives.
* Evaluated on x_i = [0.0, 1.0].
* Multiobjective fitness function. (num_objectives = num_obj)
* The optimal solutions are: sum(f) = 0.5
*/
class DTLZ1
{
public:

    explicit DTLZ1(size_t nvars = 7, size_t nobj = 3) noexcept :
        num_vars(nvars),
        num_obj(nobj) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() > num_obj);
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        std::vector<double> x1(x.begin(), x.begin() + num_obj - 1);
        std::vector<double> x2(x.begin() + num_obj - 1, x.end());

        double gm = g(x2);

        std::vector<double> fitness(num_obj, 1.0);
        
        /* F0 */
        for (const auto& val : x1)
        {
            fitness[0] *= val;
        }

        /* F1- */
        for (size_t i = 1; i < num_obj; i++)
        {
            for (size_t j = 0; j < x1.size() - i; j++)
            {
                fitness[i] *= x1[j];
            }
            fitness[i] *= (1.0 - x1[x1.size() - i]);
        }

        /* Maximization. */
        for (auto& f : fitness)
        {
            f *= -0.5 * (1.0 + gm);
        }

        return fitness;
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 7;
    size_t num_obj = 3;
    size_t var_bits = 32;
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double intval() noexcept { return ubound() - lbound(); }

private:

    static double g(const std::vector<double>& xm) noexcept
    {
        double g = double(xm.size());
        for (const auto& val : xm)
        {
            g += pow(val - 0.5, 2) - cos(20 * PI * (val - 0.5));
        }

        return 100 * g;
    }
};

/*
* Implementation of the DTLZ2 function for any number of dimensions and objectives.
* Evaluated on x_i = [0.0, 1.0].
* Multiobjective fitness function. (num_objectives = num_obj)
* The optimal solutions are: sum(f^2) = 1.0
*/
class DTLZ2
{
public:

    explicit DTLZ2(size_t nvars = 12, size_t nobj = 3) noexcept :
        num_vars(nvars),
        num_obj(nobj) {}

    /* For real chromosomes. */
    std::vector<double> operator()(const std::vector<double>& x) const
    {
        assert(x.size() > num_obj);
        assert(x.size() == num_vars);
        assert(std::all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        std::vector<double> x1(x.begin(), x.begin() + num_obj - 1);
        std::vector<double> x2(x.begin() + num_obj - 1, x.end());

        double gm = g(x2);

        std::vector<double> fitness(num_obj, 1.0);

        /* F0 */
        for (const auto& val : x1)
        {
            fitness[0] *= std::cos(val * PI / 2.0);
        }

        /* F1- */
        for (size_t i = 1; i < num_obj; i++)
        {
            for (size_t j = 0; j < x1.size() - i; j++)
            {
                fitness[i] *= std::cos(x1[j] * PI / 2.0);
            }
            fitness[i] *= std::sin(x1[x1.size() - i] * PI / 2.0);
        }

        /* Maximization. */
        for (auto& f : fitness)
        {
            f *= -(1 + gm);
        }

        return fitness;
    }

    /* For binary chromosomes. */
    std::vector<double> operator()(const std::vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        auto vars = convertToReals(x, var_bits, intval(), lbound());

        return operator()(vars);
    }

    size_t num_vars = 12;
    size_t num_obj = 3;
    size_t var_bits = 32;
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
    auto bounds() { return std::vector(num_vars, GeneBounds{ lbound(), ubound() }); }
    constexpr static double intval() noexcept { return ubound() - lbound(); }

private:

    static double g(const std::vector<double>& xm) noexcept
    {
        double g = 0.0;
        for (const auto& val : xm)
        {
            g += std::pow(val - 0.5, 2);
        }

        return g;
    }
};

/* Permutation fitness functions. */

/* Traveling salesman problem. */
template<size_t N>
class TSP
{
private:
    std::array<std::array<double, N>, N> dmat;
public:
    using Coords = std::array<double, 2>;

    TSP(const std::array<Coords, N>& cities) noexcept
    {
        for (size_t i = 0; i < dmat.size(); i++)
            for (size_t j = 0; j < dmat.size(); j++)
                dmat[i][j] = std::sqrt(std::pow(cities[i][0] - cities[j][0], 2) + std::pow(cities[i][1] - cities[j][1], 2));
    }

    std::vector<double> operator()(const std::vector<size_t>& x) const
    {
        assert(x.size() == N);

        double distance = 0.0;
        for (size_t i = 0; i < x.size() - 1; i++)
        {       
            distance += dmat[x[i]][x[i + 1]];
        }
        distance += dmat[x.front()][x.back()];

        return { -distance };
    }

    constexpr static size_t num_vars() noexcept { return N; }
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double optimal_x() noexcept { return std::numeric_limits<double>::quiet_NaN(); }
    constexpr static double optimal_value() noexcept
    {
        switch (N)
        {
            case 52:
                return -7542.0;
            case 76:
                return -108159.0;
            case 124:
                return -59030.0;
            case 152:
                return -73682.0;
            case 226:
                return -80369.0;
            case 299:
                return -48191.0;
            case 439:
                return -107217.0;
            default:
                return 0.0;
        }
    }
};

/* Integer fitness functions. */

/*
* Fitness function for testing the integerGA.
* The goal is to match a target string with the GA.
*/
class MatchString
{
public:

    explicit MatchString(std::string target) :
        target_(target),
        num_vars_(target.size()) {}

    std::vector<double> operator()(const std::vector<size_t>& x) const
    {
        assert(x.size() == num_vars());

        double fitness = 0.0;
        for (size_t i = 0; i < x.size(); i++)
        {
            fitness += double(char(x[i] + 32) == target_[i]);
        }

        return { fitness };
    }

    size_t num_vars() const noexcept { return num_vars_; }
    constexpr static size_t num_obj() noexcept { return 1; }
    double optimal_value() const noexcept { return double(num_vars_); }
    std::string optimal_x() const noexcept { return target_; }

private:
    std::string target_;
    size_t num_vars_;
};

#endif // !FITNESS_FUNCTIONS_HPP