/*
* Fitness functions for testing the genetic algorithms.
* Includes some functions for both single- and multi-objective algorithms.
* All objective functions are to be maximized using the GAs.
* Note: All of the functions are modified for maximization where needed, and return vector<double>.
*/

#ifndef FITNESS_FUNCTIONS_H
#define FITNESS_FUNCTIONS_H

#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <cassert>

using namespace std;

constexpr double PI = 3.14159265358979323846;

/* Utility functions. */

/* Convert binary vector to real values in [0.0, upper_limit] for the binary fitness functions. */
vector<double> convertToReals(const vector<char>& binv, size_t bits_per_var, double upper_limit)
{
    size_t var_count = size_t(double(binv.size()) / bits_per_var);
    
    vector<double> vars(var_count);
    for (size_t i = 0; i < var_count; i++)
    {
        /* Convert to decimal. */
        size_t k = accumulate(binv.begin() + i * bits_per_var, binv.begin() + (i + 1) * bits_per_var, size_t{ 0 },
        [](size_t x, size_t y)
        {
            return (x << 1) + y;
        });

        /* Norm decimal. */
        vars[i] = double(k) / ((1ULL << bits_per_var) - 1.0);

        /* Scale. */
        vars[i] *= upper_limit;
    }

    return vars;
}


/* Single-objective fitness functions. */

/*
* Implementation of the Rastrigin function for any number of dimensions.
* Evaluated on x_i = [-5.12, 5.12].
* The global optimum of the function is f(x) = 0, at x = (0, 0, ... , 0).
*/
class Rastrigin
{
public:

    explicit Rastrigin(size_t num_vars = 10) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double fx = 10.0 * x.size();
        for (size_t i = 0; i < x.size(); i++)
        {
            fx += pow(x[i], 2) - 10.0 * cos(2 * PI * x[i]);
        }

        return { -fx };	/* For maximization. */
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -5.12; }
    constexpr static double ubound() noexcept { return 5.12; }
    constexpr static double optimal_value() noexcept { return 0.0; }
    constexpr static double optimal_x() noexcept { return 0.0; }
    constexpr static double intval() noexcept { return ubound() - lbound(); }
};

/*
* Implementation of the Rosenbrock function for any dimensions.
* Evaluated on x_i = [-2.048, 2.048].
* The global optimum of the function is f(x) = 0, at x = (1, 1, ... , 1).
*/
class Rosenbrock
{
public:

    Rosenbrock(size_t num_vars = 3) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double fx = 0.0;
        for (size_t i = 0; i < x.size() - 1; i++)
        {
            fx += 100.0 * pow((x[i + 1] - pow(x[i], 2)), 2) + pow(1 - x[i], 2);
        }

        return { -fx }; /* For maximization. */
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars =  3;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -2.048; }
    constexpr static double ubound() noexcept { return 2.048; }
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

    explicit Schwefel(size_t num_vars = 10) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double fx = 418.9829 * x.size();
        for (size_t i = 0; i < x.size(); i++)
        {
            fx -= x[i] * sin(sqrt(abs(x[i])));
        }

        return { -fx }; /* For maximization. */
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -500.0; }
    constexpr static double ubound() noexcept { return 500.0; }
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

    explicit Griewank(size_t num_vars = 10) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double fx = 1.0;
        double fminus = 1.0;
        for (size_t i = 0; i < x.size(); i++)
        {
            fx += x[i] * x[i] / 4000.0;
            fminus *= cos(x[i] / sqrt(i + 1.0));
        }
        fx -= fminus;

        return { -fx }; /* For maximization. */
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -600.0; }
    constexpr static double ubound() noexcept { return 600.0; }
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

    explicit Ackley(size_t num_vars = 10) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = 0.0;
        double f2 = 0.0;
        for (size_t i = 0; i < x.size(); i++)
        {
            f1 += pow(x[i], 2);
            f2 += cos(2.0 * PI * x[i]);
        }
        f1 = exp(-0.2 * sqrt(f1 / num_vars));
        f2 = exp(f2 / num_vars);
        
        double fx = -20.0 * f1 - f2 + 20.0 + exp(1.0);

        return { -fx }; /* For maximization. */
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 1; }
    constexpr static double lbound() noexcept { return -32.768; }
    constexpr static double ubound() noexcept { return 32.768; }
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

    explicit KUR(size_t num_vars = 3) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() > 1);
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = 0.0;
        for (size_t i = 0; i < x.size() - 1; i++)
        {
            f1 += -10.0 * exp(-0.2 * sqrt(pow(x[i], 2) + pow(x[i + 1], 2)));
        }
        double f2 = 0.0;
        for (size_t i = 0; i < x.size(); i++)
        {
            f2 += pow(abs(x[i]), 0.8) + 5 * sin(pow(x[i], 3));
        }

        return { -f1, -f2 };
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars = 3;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 2; }
    constexpr static double lbound() noexcept { return -5.0; }
    constexpr static double ubound() noexcept { return 5.0; }
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

    explicit ZDT2(size_t num_vars = 30) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() > 1);
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = x[0];
        double g = accumulate(x.begin() + 1, x.end(), 0.0);
        g = 9 * g / (x.size() - 1.0) + 1;
        double f2 = g * (1.0 - pow(x[0] / g, 2));

        return { -f1, -f2 };
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars = 30;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 2; }
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
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

    explicit ZDT3(size_t num_vars = 30) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() > 1);
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = x[0];
        double g = accumulate(x.begin() + 1, x.end(), 0.0);
        g = 1 + 9 * g / (x.size() - 1.0);
        double f2 = g * (1 - sqrt(x[0] / g) - (x[0] / g) * sin(10 * PI * x[0]));

        return { -f1, -f2 };
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars = 30;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 2; }
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
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

    explicit ZDT6(size_t num_vars = 10) : num_vars(num_vars) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() > 1);
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        double f1 = 1 - exp(-4 * x[0]) * pow(sin(6 * PI * x[0]), 6);
        double g = accumulate(x.begin() + 1, x.end(), 0.0);
        g = 9.0 * pow(g / (x.size() - 1.0), 0.25) + 1;
        double f2 = g * (1 - pow(f1 / g, 2));

        return { -f1, -f2 };
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars = 10;
    size_t var_bits = 32;
    constexpr static size_t num_obj() noexcept { return 2; }
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
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

    explicit DTLZ1(size_t num_vars = 7, size_t num_obj = 3) : num_vars(num_vars), num_obj(num_obj) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() > num_obj);
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        vector<double> x1(x.begin(), x.begin() + num_obj - 1);
        vector<double> x2(x.begin() + num_obj - 1, x.end());

        double gm = g(x2);

        vector<double> fitness(num_obj, 1.0);
        
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
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars;
    size_t num_obj;
    size_t var_bits = 32;
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
    constexpr static double intval() noexcept { return ubound() - lbound(); }

private:

    static double g(const vector<double>& xm)
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

    explicit DTLZ2(size_t num_vars = 12, size_t num_obj = 3) : num_vars(num_vars), num_obj(num_obj) {}

    /* For real chromosomes. */
    vector<double> operator()(const vector<double>& x) const
    {
        assert(x.size() > num_obj);
        assert(x.size() == num_vars);
        assert(all_of(x.begin(), x.end(), [](double val) { return lbound() <= val && val <= ubound(); }));

        vector<double> x1(x.begin(), x.begin() + num_obj - 1);
        vector<double> x2(x.begin() + num_obj - 1, x.end());

        double gm = g(x2);

        vector<double> fitness(num_obj, 1.0);

        /* F0 */
        for (const auto& val : x1)
        {
            fitness[0] *= cos(val * PI / 2.0);
        }

        /* F1- */
        for (size_t i = 1; i < num_obj; i++)
        {
            for (size_t j = 0; j < x1.size() - i; j++)
            {
                fitness[i] *= cos(x1[j] * PI / 2.0);
            }
            fitness[i] *= sin(x1[x1.size() - i] * PI / 2.0);
        }

        /* Maximization. */
        for (auto& f : fitness)
        {
            f *= -(1 + gm);
        }

        return fitness;
    }

    /* For binary chromosomes. */
    vector<double> operator()(const vector<char>& x) const
    {
        assert(x.size() == num_vars * var_bits);

        /* Convert binary chromosome to real values. */
        vector<double> vars = convertToReals(x, var_bits, intval());
        for_each(vars.begin(), vars.end(), [](double& var) { var += lbound(); });

        return operator()(vars);
    }

    size_t num_vars;
    size_t num_obj;
    size_t var_bits = 32;
    constexpr static double lbound() noexcept { return 0.0; }
    constexpr static double ubound() noexcept { return 1.0; }
    constexpr static double intval() noexcept { return ubound() - lbound(); }

private:

    static double g(const vector<double>& xm)
    {
        double g = 0.0;
        for (const auto& val : xm)
        {
            g += pow(val - 0.5, 2);
        }

        return g;
    }
};

/* Permutation fitness functions. */

/*
* Traveling salesman problem. The node coordinates are read from files.
*/
class TSP
{
private:

    size_t num_vars_ = 0;
    vector<pair<double, double>> coords;

public:

    TSP(string fname)
    {
        ifstream file(fname);

        string line;
        while (getline(file, line))
        {
            stringstream linestream(line);
            int idx;
            double c1, c2;

            linestream >> idx >> c1 >> c2;
            coords.emplace_back(c1, c2);
        }
        file.close();

        num_vars_ = coords.size();
    }

    vector<double> operator()(const vector<size_t>& x) const
    {
        assert(x.size() == num_vars());

        double distance = 0.0;
        for (size_t i = 0; i < x.size() - 1; i++)
        {
            pair<double, double> node1 = coords[x[i]];
            pair<double, double> node2 = coords[x[i + 1]];
            
            distance += sqrt(pow(node1.first - node2.first, 2) + pow(node1.second - node2.second, 2));
        }
        pair<double, double> node1 = coords[x.front()];
        pair<double, double> node2 = coords[x.back()];

        distance += sqrt(pow(node1.first - node2.first, 2) + pow(node1.second - node2.second, 2));

        return { -distance };
    }

    size_t num_vars() const noexcept { return num_vars_; }
    constexpr static size_t num_obj() noexcept { return 1; }
    double optimal_value() const noexcept
    {
        switch (num_vars_)
        {
            case 52:
                return -7542.0;
            case 124:
                return -59030.0;
            case 226:
                return -80369.0;
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

    explicit MatchString(string target) : target_(target), num_vars_(target.size()) {}

    vector<double> operator()(const vector<size_t>& x) const
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
    string optimal_x() const noexcept { return target_; }

private:
    string target_;
    size_t num_vars_;
};

#endif // !FITNESS_FUNCTIONS_H