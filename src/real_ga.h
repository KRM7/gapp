/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this softwareand associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright noticeand this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

/**
* This file contains the real encoded genetic algorithm class.
*
* @file real_ga.h
*/

#ifndef GA_RCGA_H
#define GA_RCGA_H

#include <vector>
#include <utility>

#include "base_ga.h"

namespace genetic_algorithm
{
    /**
    * Standard genetic algorithm that uses real encoding. \n
    * Each gene of the chromosomes is a real value.
    */
    class RCGA : public GA<double>
    {
    public:

        /**
        * For the gene boundaries. \n
        * For example: { {gene1_min, gene1_max}, {gene2_min, gene2_max}, ... }
        */
        using limits_t = std::vector<std::pair<double, double>>;

        /**
        * Possible crossover operators that can be used in the RCGA. \n
        * Includes commonly used crossover operators in real-coded genetic algorithms, but a custom function
        * can also be used to perform the crossovers with the custom option. \n
        * Set the crossover method used in the algorithm with @ref crossover_method.
        */
        enum class CrossoverMethod
        {
            arithmetic,			/**< Arithmetic crossover operator.	Uses no parameters. */
            blx_a,				/**< BLX-alpha (blend) crossover operator. @see blx_crossover_param */
            simulated_binary,	/**< Simulated binary crossover (SBX) operator. @see sim_binary_crossover_param */
            wright,				/**< Wright heuristic crossover (HX) operator. Uses no parameters. */
            custom				/**< Custom crossover operator defined by the user. @see setCrossoverFunction */
        };

        /**
        * Possible mutation operators that can be used in the RCGA. \n
        * Includes commonly used mutation operators in real-coded genetic algorithms, but a custom function
        * can also be used to perform the mutations with the custom option. \n
        * Set the mutation method used in the algorithm with @ref mutation_method.
        */
        enum class MutationMethod
        {
            random,			/**<  Random (uniform) mutation operator. Uses no parameters. */
            polynomial,		/**<  Polynomial mutation operator. @see polynomial_mutation_param */
            nonuniform,		/**<  Non-uniform mutation operator. @see nonuniform_mutation_param */
            boundary,		/**<  Boundary mutation operator. Uses no parameters. */
            gauss,			/**<  Gauss mutation operator. @see gauss_mutation_param */
            custom			/**<  Custom mutation operator defined by the user. Uses the @ref customMutate to perform the mutations. */
        };

        /**
        * Basic constructor for the RCGA.
        *
        * @param chrom_len The number of real genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries of the real coded genes (their min and max values).
        */
        RCGA(size_t chrom_len, fitnessFunction_t fitnessFunction, limits_t bounds);

        /**
        * Sets the crossover function used in the algorithm to @f.
        * @see CrossoverMethod
        *
        * @param method The crossover function to use.
        */
        void crossover_method(crossoverFunction_t f);

        /**
        * Sets the crossover method used in the algorithm to @p method.
        * @see CrossoverMethod
        *
        * @param method The crossover method to use.
        */
        void crossover_method(CrossoverMethod method);
        [[nodiscard]] CrossoverMethod crossover_method() const;

        /**
        * Sets the mutation function used in the algorithm to @f.
        * @see MutationMethod
        *
        * @param method The mutation function to use.
        */
        void mutation_method(mutationFunction_t f);

        /**
        * Sets the mutation method used in the algorithm to @p method.
        * @see MutationMethod
        *
        * @param method The mutation method to use.
        */
        void mutation_method(MutationMethod method);
        [[nodiscard]] MutationMethod mutation_method() const;

        /**
        * Sets the boundaries of the real-coded genes. \n
        * Each element must contain the lower and upper bounds of the corresponding real gene. (The min and max values of the gene.) \n
        * The number of elements must be the same as the length of the chromosomes, and the lower bounds must not be higher than the upper bounds. \n
        * Eg. in the case of chromosomes of 2 values where both genes must be between -1 and 1: \n
        * limits = { {-1.0, 1.0}, {-1.0, 1.0} }
        *
        * @param limits The boundaries of the real genes.
        */
        void limits(limits_t limits);
        [[nodiscard]] limits_t limits() const;

        /**
        * Sets the alpha parameter of the BLX-alpha crossover method. \n
        * The parameter controls the length of the interval the children's genes are randomly chosen from, larger alpha -> larger interval. \n
        * Must be nonnegative, with the ideal being around 0.5.
        * @see crossover_method @see CrossoverMethod
        *
        * @param alpha The alpha parameter of the BLX-alpha crossover.
        */
        void blx_crossover_param(double alpha);
        [[nodiscard]] double blx_crossover_param() const;

        /**
        * Sets the shape parameter (eta) of the simulated binary crossover method. \n
        * Must be nonnegative, typical values are 1-5. \n
        * Larger values will lead to the children's genes to be closer to the parents' genes.
        * @see crossover_method @see CrossoverMethod
        *
        * @param eta The parameter of the simulated binary crossover.
        */
        void sim_binary_crossover_param(double eta);
        [[nodiscard]] double sim_binary_crossover_param() const;

        /**
        * Sets the time parameter of the non-uniform mutation operator. \n
        * The mutated genes will be closer to the original values as the generations advance. Larger parameter values
        * will result in this process happening faster. \n
        * A value of 0 means the behaviour of the mutation does not change over time. \n
        * Must be nonnegative.
        * @see mutation_method @see MutationMethod
        *
        * @param b The parameter of the non-uniform mutation.
        */
        void nonuniform_mutation_param(double b);
        [[nodiscard]] double nonuniform_mutation_param() const;

        /**
        * Sets the parameter of the polynomial mutation operator. \n
        * Must be nonnegative, typically in the range [20, 100]. \n
        * The mutated genes will be closer to their original values with higher parameter values.
        * @see mutation_method @see MutationMethod
        *
        * @param eta The parameter of the polynomial mutation.
        */
        void polynomial_mutation_param(double eta);
        [[nodiscard]] double polynomial_mutation_param() const;

        /**
        * Sets the parameter of the gauss mutation operator. \n
        * Must be larger than 0. \n
        * Controls how many standard deviations the gene's allowed range should take up. Larger values will result in mutated genes
        * closer to the original genes.
        * @see mutation_method @see MutationMethod
        *
        * @param sigmas The parameter of the gauss mutation.
        */
        void gauss_mutation_param(double sigmas);
        [[nodiscard]] double gauss_mutation_param() const;

    private:

        limits_t limits_;

        CrossoverMethod crossover_method_ = CrossoverMethod::blx_a;
        double blx_crossover_param_ = 0.5;
        double sim_binary_crossover_param_ = 4.0;

        MutationMethod mutation_method_ = MutationMethod::random;
        double nonuniform_mutation_param_ = 2.0;
        double polynomial_mutation_param_ = 40.0;
        double gauss_mutation_param_ = 6.0;

        Candidate generateCandidate() const override;
        CandidatePair crossover(const Candidate& p1, const Candidate& p2) const override;
        void mutate(Candidate& child) const override;

        static CandidatePair arithmeticCrossover(const Candidate& parent1, const Candidate& parent2, double pc);
        static CandidatePair blxAlphaCrossover(const Candidate& parent1, const Candidate& parent2, double pc, double alpha, const limits_t& bounds);
        static CandidatePair simulatedBinaryCrossover(const Candidate& parent1, const Candidate& parent2, double pc, double b, const limits_t& bounds);
        static CandidatePair wrightCrossover(const Candidate& parent1, const Candidate& parent2, double pc, const limits_t& bounds);

        static void randomMutate(Candidate& child, double pm, const limits_t& bounds);
        static void nonuniformMutate(Candidate& child, double pm, size_t time, size_t time_max, double b, const limits_t& bounds);
        static void polynomialMutate(Candidate& child, double pm, double eta, const limits_t& bounds);
        static void boundaryMutate(Candidate& child, double pm, const limits_t& bounds);
        static void gaussMutate(Candidate& child, double pm, double scale, const limits_t& bounds);
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cassert>

#include "rng.h"

namespace genetic_algorithm
{
    inline RCGA::RCGA(size_t chrom_len, fitnessFunction_t fitnessFunction, limits_t bounds)
        : GA(chrom_len, fitnessFunction), limits_(bounds)
    {
        if (bounds.size() != chrom_len)
        {
            throw std::invalid_argument("The size of the bounds must be the same as the number of genes.");
        }
        if (std::any_of(bounds.begin(), bounds.end(), [](std::pair<double, double> b) { return b.first > b.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

    }

    template<>
    inline bool RCGA::Candidate::operator==(const Candidate& rhs) const
    {
        return std::equal(this->chromosome.begin(), this->chromosome.end(), rhs.chromosome.begin(),
        [](double lhs, double rhs)
        {
            return std::abs(lhs - rhs) <= std::numeric_limits<double>::epsilon() * std::max(std::abs(lhs), std::abs(rhs));
        });
    }

    template<>
    inline bool RCGA::Candidate::operator!=(const Candidate& rhs) const
    {
        return !std::equal(this->chromosome.begin(), this->chromosome.end(), rhs.chromosome.begin(),
        [](double lhs, double rhs)
        {
            return std::abs(lhs - rhs) <= std::numeric_limits<double>::epsilon() * std::max(std::abs(lhs), std::abs(rhs));
        });
    }


    inline void RCGA::crossover_method(crossoverFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        crossover_method_ = CrossoverMethod::custom;
        customCrossover = f;
    }

    inline void RCGA::crossover_method(CrossoverMethod method)
    {
        if (static_cast<size_t>(method) > 4) throw std::invalid_argument("Invalid crossover method selected.");

        crossover_method_ = method;
    }

    inline RCGA::CrossoverMethod RCGA::crossover_method() const
    {
        return crossover_method_;
    }

    inline void RCGA::mutation_method(mutationFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        mutation_method_ = MutationMethod::custom;
        customMutate = f;
    }

    inline void RCGA::mutation_method(MutationMethod method)
    {
        if (static_cast<size_t>(method) > 5) throw std::invalid_argument("Invalid mutation method selected.");

        mutation_method_ = method;
    }

    inline RCGA::MutationMethod RCGA::mutation_method() const
    {
        return mutation_method_;
    }

    inline void RCGA::limits(limits_t limits)
    {
        if (limits.size() != chrom_len_)
        {
            throw std::invalid_argument("The number of limits must be equal to the chromosome length.");
        }
        if (std::any_of(limits.begin(), limits.end(), [](std::pair<double, double> b) {return b.first > b.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

        limits_ = limits;
    }

    inline RCGA::limits_t RCGA::limits() const
    {
        return limits_;
    }

    inline void RCGA::blx_crossover_param(double alpha)
    {
        if (!(0.0 <= alpha && alpha <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Alpha must be a nonnegative, finite value.");
        }

        blx_crossover_param_ = alpha;
    }

    inline double RCGA::blx_crossover_param() const
    {
        return blx_crossover_param_;
    }

    inline void RCGA::sim_binary_crossover_param(double eta)
    {
        if (!(0.0 <= eta && eta <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Eta must be a nonnegative, finite value.");
        }

        sim_binary_crossover_param_ = eta;
    }

    inline double RCGA::sim_binary_crossover_param() const
    {
        return sim_binary_crossover_param_;
    }

    inline void RCGA::nonuniform_mutation_param(double b)
    {
        if (!(0.0 <= b && b <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("The parameter b must be a nonnegative, finite value.");
        }

        nonuniform_mutation_param_ = b;
    }

    inline double RCGA::nonuniform_mutation_param() const
    {
        return nonuniform_mutation_param_;
    }

    inline void RCGA::polynomial_mutation_param(double eta)
    {
        if (!(0.0 <= eta && eta <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Eta must be a nonnegative, finite value.");
        }

        polynomial_mutation_param_ = eta;
    }

    inline double RCGA::polynomial_mutation_param() const
    {
        return polynomial_mutation_param_;
    }

    inline void RCGA::gauss_mutation_param(double sigmas)
    {
        if (!(0.0 < sigmas && sigmas <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Eta must be a nonnegative, finite value.");
        }

        gauss_mutation_param_ = sigmas;
    }

    inline double RCGA::gauss_mutation_param() const
    {
        return gauss_mutation_param_;
    }


    inline RCGA::Candidate RCGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);
        assert(chrom_len_ == limits_.size());

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(rng::randomReal(limits_[i].first, limits_[i].second));
        }

        return sol;
    }

    inline RCGA::CandidatePair RCGA::crossover(const Candidate& p1, const Candidate& p2) const
    {
        switch (crossover_method_)
        {
            case CrossoverMethod::arithmetic:
                return arithmeticCrossover(p1, p2, crossover_rate_);
            case CrossoverMethod::blx_a:
                return blxAlphaCrossover(p1, p2, crossover_rate_, blx_crossover_param_, limits_);
            case CrossoverMethod::simulated_binary:
                return simulatedBinaryCrossover(p1, p2, crossover_rate_, sim_binary_crossover_param_, limits_);
            case CrossoverMethod::wright:
                return wrightCrossover(p1, p2, crossover_rate_, limits_);
            case CrossoverMethod::custom:
                return customCrossover(p1, p2, crossover_rate_);
            default:
                assert(false);	/* Invalid crossover method. Shouldn't get here. */
                std::abort();
        }
    }

    inline void RCGA::mutate(Candidate& child) const
    {
        switch (mutation_method_)
        {
            case MutationMethod::random:
                randomMutate(child, mutation_rate_, limits_);
                break;
            case MutationMethod::nonuniform:
                nonuniformMutate(child, mutation_rate_, generation_cntr_, max_gen_, nonuniform_mutation_param_, limits_);
                break;
            case MutationMethod::polynomial:
                polynomialMutate(child, mutation_rate_, polynomial_mutation_param_, limits_);
                break;
            case MutationMethod::boundary:
                boundaryMutate(child, mutation_rate_, limits_);
                break;
            case MutationMethod::gauss:
                gaussMutate(child, mutation_rate_, gauss_mutation_param_, limits_);
                break;
            case MutationMethod::custom:
                customMutate(child, mutation_rate_);
                break;
            default:
                assert(false);	/* Invalid mutation method. Shouldnt get here. */
                std::abort();
        }
    }

    inline RCGA::CandidatePair RCGA::arithmeticCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
    {
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(0.0 <= pc && pc <= 1.0);

        Candidate child1(parent1), child2(parent2);

        /* Perform crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            double alpha = rng::randomReal();
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                child1.chromosome[i] = alpha * parent1.chromosome[i] + (1.0 - alpha) * parent2.chromosome[i];
                child2.chromosome[i] = (1.0 - alpha) * parent1.chromosome[i] + alpha * parent2.chromosome[i];
            }
            child1.is_evaluated = false;
            child2.is_evaluated = false;
        }

        return std::make_pair(child1, child2);
    }

    inline RCGA::CandidatePair RCGA::blxAlphaCrossover(const Candidate& parent1, const Candidate& parent2, double pc, double alpha, const limits_t& bounds)
    {
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(parent1.chromosome.size() == bounds.size());
        assert(0.0 <= pc && pc <= 1.0);
        assert(alpha >= 0.0);

        Candidate child1(parent1), child2(parent2);

        /* Perform crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* Calc interval to generate the childrens genes on. */
                auto [range_min, range_max] = std::minmax(parent1.chromosome[i], parent2.chromosome[i]);
                double range_ext = alpha * (range_max - range_min);
                /* Generate genes from an uniform distribution on the interval. */
                child1.chromosome[i] = rng::randomReal(range_min - range_ext, range_max + range_ext);
                child2.chromosome[i] = rng::randomReal(range_min - range_ext, range_max + range_ext);
                /* The generated genes might be outside the allowed interval. */
                child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].first, bounds[i].second);
                child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].first, bounds[i].second);
            }
            child1.is_evaluated = false;
            child2.is_evaluated = false;
        }

        return std::make_pair(child1, child2);
    }

    inline RCGA::CandidatePair RCGA::simulatedBinaryCrossover(const Candidate& parent1, const Candidate& parent2, double pc, double b, const limits_t& bounds)
    {
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(parent1.chromosome.size() == bounds.size());
        assert(0.0 <= pc && pc <= 1.0);
        assert(b > 0.0);

        Candidate child1(parent1), child2(parent2);

        /* Perform crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            /* Generate beta from the distribution. */
            double u = rng::randomReal();
            double beta = (u <= 0.5) ? std::pow(2.0 * u, 1.0 / (b + 1.0)) : std::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (b + 1.0));

            /* Perform crossover. */
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                child1.chromosome[i] = 0.5 * ((1 - beta) * parent1.chromosome[i] + (1 + beta) * parent2.chromosome[i]);
                child2.chromosome[i] = 0.5 * ((1 + beta) * parent1.chromosome[i] + (1 - beta) * parent2.chromosome[i]);
                /* The childrens genes may be outside the allowed interval. */
                child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].first, bounds[i].second);
                child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].first, bounds[i].second);
            }
            child1.is_evaluated = false;
            child2.is_evaluated = false;
        }

        return std::make_pair(child1, child2);
    }

    inline RCGA::CandidatePair RCGA::wrightCrossover(const Candidate& parent1, const Candidate& parent2, double pc, const limits_t& bounds)
    {
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(parent1.chromosome.size() == bounds.size());
        assert(0.0 <= pc && pc <= 1.0);

        Candidate child1(parent1), child2(parent2);

        /* Perform crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            /* p1 is always the better parent. */
            const Candidate* p1 = detail::paretoCompare(parent1.fitness, parent2.fitness) ? &parent2 : &parent1;
            const Candidate* p2 = detail::paretoCompare(parent1.fitness, parent2.fitness) ? &parent1 : &parent2;
            /* Get random weights. */
            double w1 = rng::randomReal();
            double w2 = rng::randomReal();
            /* Perform crossover. */
            for (size_t i = 0; i < p1->chromosome.size(); i++)
            {
                child1.chromosome[i] = w1 * (p1->chromosome[i] - p2->chromosome[i]) + p1->chromosome[i];
                child2.chromosome[i] = w2 * (p1->chromosome[i] - p2->chromosome[i]) + p1->chromosome[i];
                /* The generated childrens genes may be outside the allowed intervals. */
                child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].first, bounds[i].second);
                child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].first, bounds[i].second);
            }
            child1.is_evaluated = false;
            child2.is_evaluated = false;
        }

        return std::make_pair(child1, child2);
    }

    inline void RCGA::randomMutate(Candidate& child, double pm, const limits_t& bounds)
    {
        assert(0.0 <= pm && pm <= 1.0);
        assert(child.chromosome.size() == bounds.size());

        for (size_t i = 0; i < child.chromosome.size(); i++)
        {
            /* Mutate the gene with pm probability. */
            if (rng::randomReal() <= pm)
            {
                child.chromosome[i] = rng::randomReal(bounds[i].first, bounds[i].second);
                child.is_evaluated = false;
            }
        }
    }

    inline void RCGA::nonuniformMutate(Candidate& child, double pm, size_t time, size_t time_max, double b, const limits_t& bounds)
    {
        assert(0.0 <= pm && pm <= 1.0);
        assert(child.chromosome.size() == bounds.size());
        assert(b >= 0.0);

        for (size_t i = 0; i < child.chromosome.size(); i++)
        {
            /* Perform mutation on the gene with pm probability. */
            if (rng::randomReal() <= pm)
            {
                double interval = bounds[i].second - bounds[i].first;
                double r = rng::randomReal();
                double sign = rng::randomBool() ? 1.0 : -1.0;

                child.chromosome[i] += sign * interval * (1.0 - std::pow(r, std::pow(1.0 - double(time) / time_max, b)));
                child.is_evaluated = false;

                /* The mutated gene might be outside the allowed range. */
                child.chromosome[i] = std::clamp(child.chromosome[i], bounds[i].first, bounds[i].second);
            }
        }
    }

    inline void RCGA::polynomialMutate(Candidate& child, double pm, double eta, const limits_t& bounds)
    {
        assert(0.0 <= pm && pm <= 1.0);
        assert(child.chromosome.size() == bounds.size());
        assert(eta >= 0.0);

        for (size_t i = 0; i < child.chromosome.size(); i++)
        {
            /* Perform mutation on the gene with pm probability. */
            if (rng::randomReal() <= pm)
            {
                double u = rng::randomReal();
                if (u <= 0.5)
                {
                    double delta = std::pow(2.0 * u, 1.0 / (1.0 + eta)) - 1.0;
                    child.chromosome[i] += delta * (child.chromosome[i] - bounds[i].first);
                }
                else
                {
                    double delta = 1.0 - std::pow(2.0 - 2.0 * u, 1.0 / (1.0 + eta));
                    child.chromosome[i] += delta * (bounds[i].second - child.chromosome[i]);
                }
                child.is_evaluated = false;
                /* The mutated gene will always be in the allowed range. */
            }
        }
    }

    inline void RCGA::boundaryMutate(Candidate& child, double pm, const limits_t& bounds)
    {
        assert(0.0 <= pm && pm <= 1.0);
        assert(child.chromosome.size() == bounds.size());

        for (size_t i = 0; i < child.chromosome.size(); i++)
        {
            /* Perform mutation on the gene with pm probability. */
            if (rng::randomReal() <= pm)
            {
                child.chromosome[i] = rng::randomBool() ? bounds[i].first : bounds[i].second;
                child.is_evaluated = false;
            }
        }
    }

    inline void RCGA::gaussMutate(Candidate& child, double pm, double scale, const limits_t& bounds)
    {
        assert(0.0 <= pm && pm <= 1.0);
        assert(child.chromosome.size() == bounds.size());
        assert(scale > 0.0);

        for (size_t i = 0; i < child.chromosome.size(); i++)
        {
            /* Perform mutation on the gene with pm probability. */
            if (rng::randomReal() <= pm)
            {
                double SD = (bounds[i].second - bounds[i].first) / scale;
                child.chromosome[i] += rng::randomNormal(0.0, SD);
                child.is_evaluated = false;
                /* The mutated gene might be outside the allowed range. */
                child.chromosome[i] = std::clamp(child.chromosome[i], bounds[i].first, bounds[i].second);
            }
        }
    }

} // namespace genetic_algorithm

#endif // !GA_RCGA_H