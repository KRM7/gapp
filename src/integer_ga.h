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
* This file contains the integer coded genetic algorithm class.
*
* @file integer_ga.h
*/

#ifndef GA_INTEGER_GA_H
#define GA_INTEGER_GA_H

#include <cstddef>

#include "base_ga.h"

namespace genetic_algorithm
{
    /**
    * Integer coded GA. \n
    * Same as @ref BinaryGA, but the genes of the chromosomes can be any integer on [0, base], not just 0 or 1. \n
    * It also uses a slightly different mutation function with swaps and inversions.
    */
    class IntegerGA : public GA<size_t>
    {
    public:

        /**
        * Possible crossover operators that can be used in the IntegerGA. \n
        * These crossover method are the same as the ones used in the binary coded algorithm. \n
        * Set the crossover method used in the algorithm with @ref crossover_method. \n
        * The function used for the crossovers with the custom method can be set with @ref setCrossoverFunction.
        */
        enum class CrossoverMethod
        {
            single_point,	/**< Single-point crossover operator. */
            two_point,		/**< Two-point crossover operator. */
            n_point,		/**< General n-point crossover operator. Set the number of crossover points with @ref num_crossover_points */
            uniform,		/**< Uniform crossover operator. */
            custom			/**< Custom crossover operator defined by the user. @see setCrossoverMethod */
        };

        /**
        * Possible mutation operators that can be used in the IntegerGA. \n
        * Same operators as in the binary coded algorithm, with some small changes. \n
        * Set the mutation method used in the algorithm with @ref mutation_method. \n
        * The function used for the mutations with the custom method can be set with @ref setMutationFunction.
        */
        enum class MutationMethod
        {
            standard,	/**< Standard mutation operator used in the @ref BinaryGA with swap and inversion added. @see swap_rate @see inversion_rate */
            custom		/**< Custom mutation operator defined by the user. @see setMutationFunction */
        };

        /**
        * Basic contructor for the IntegerGA.
        *
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function used in the algorithm.
        * @param base The number of values a gene can take. Must be > 1. If 2, same as the @ref BinaryGA.
        */
        IntegerGA(size_t chrom_len, fitnessFunction_t fitnessFunction, size_t base);

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
        * Sets the number of crossover points used in the crossovers to @p n if the n_point crossover method selected. \n
        * The number of crossover points must be at least 1.
        * @see crossover_method @see CrossoverMethod
        *
        * @param n The number of crossover points.
        */
        void num_crossover_points(size_t n);
        [[nodiscard]] size_t num_crossover_points() const;

        /**
        * Sets the number of values a gene can take to @p base. \n
        * The value of the base must be at least 2, and the GA is essentially the same as the
        * BinaryGA if the base is set to 2.
        *
        * @param base The number of values a gene can be.
        */
        void base(size_t base);
        [[nodiscard]] size_t base() const;

        /**
        * Sets the probability of a single swap occuring during the mutation of a Candidate to @p ps. \n
        * The value of ps must be on the closed interval [0.0, 1.0].
        *
        * @param ps The probability of swap during mutation.
        */
        void swap_rate(double ps);
        [[nodiscard]] double swap_rate() const;

        /**
        * Sets the probability of inversion during the mutation of a Candidate to @p pi. \n
        * The value of pi must be on the closed interval [0.0, 1.0].
        *
        * @param pi The probability of inversion during mutation.
        */
        void inversion_rate(double pi);
        [[nodiscard]] double inversion_rate() const;

    private:

        CrossoverMethod crossover_method_ = CrossoverMethod::single_point;
        MutationMethod mutation_method_ = MutationMethod::standard;
        size_t num_crossover_points_ = 3;
        size_t base_ = 4;
        double swap_rate_ = 0.1;
        double inversion_rate_ = 0.1;

        Candidate generateCandidate() const override;
        CandidatePair crossover(const Candidate& parent1, const Candidate& parent2) const override;
        void mutate(Candidate& child) const override;

        static CandidatePair nPointCrossover(const Candidate& parent1, const Candidate& parent2, double pc, size_t n);
        static CandidatePair uniformCrossover(const Candidate& parent1, const Candidate& parent2, double pc);

        static void standardMutate(Candidate& child, double pm, double ps, double pi, size_t base_);
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <utility>
#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "rng.h"

namespace genetic_algorithm
{
    inline IntegerGA::IntegerGA(size_t chrom_len, fitnessFunction_t fitnessFunction, size_t base)
        : GA(chrom_len, fitnessFunction), base_(base)
    {
        if (base < 2) throw std::invalid_argument("The base must be at least 2.");
    }


    inline void IntegerGA::crossover_method(crossoverFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        crossover_method_ = CrossoverMethod::custom;
        customCrossover = f;
    }

    inline void IntegerGA::crossover_method(CrossoverMethod method)
    {
        if (static_cast<size_t>(method) > 4) throw std::invalid_argument("Invalid crossover method selected.");

        crossover_method_ = method;
    }

    inline IntegerGA::CrossoverMethod IntegerGA::crossover_method() const
    {
        return crossover_method_;
    }

    inline void IntegerGA::mutation_method(mutationFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        mutation_method_ = MutationMethod::custom;
        customMutate = f;
    }

    inline void IntegerGA::mutation_method(MutationMethod method)
    {
        if (static_cast<size_t>(method) > 1) throw std::invalid_argument("Invalid mutation method selected.");

        mutation_method_ = method;
    }

    inline IntegerGA::MutationMethod IntegerGA::mutation_method() const
    {
        return mutation_method_;
    }

    inline void IntegerGA::num_crossover_points(size_t n)
    {
        if (n == 0) throw std::invalid_argument("The number of crossover points must be at least 1.");

        num_crossover_points_ = n;
    }

    inline size_t IntegerGA::num_crossover_points() const
    {
        return num_crossover_points_;
    }

    inline void IntegerGA::base(size_t base)
    {
        if (base < 2) throw std::invalid_argument("The base must be at least 2.");

        base_ = base;
    }

    inline size_t IntegerGA::base() const
    {
        return base_;
    }

    inline void IntegerGA::swap_rate(double ps)
    {
        if (!(0.0 <= ps && ps <= 1.0)) throw std::invalid_argument("The probability of a swap must be in [0, 1].");

        swap_rate_ = ps;
    }

    inline double IntegerGA::swap_rate() const
    {
        return swap_rate_;
    }

    inline void IntegerGA::inversion_rate(double pi)
    {
        if (!(0.0 <= pi && pi <= 1.0)) throw std::invalid_argument("The probability of inversion must be in [0, 1].");

        inversion_rate_ = pi;
    }

    inline double IntegerGA::inversion_rate() const
    {
        return inversion_rate_;
    }


    inline IntegerGA::Candidate IntegerGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);
        assert(base_ > 1);

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(rng::randomInt(size_t{ 0 }, base_ - 1));
        }

        return sol;
    }

    inline IntegerGA::CandidatePair IntegerGA::crossover(const Candidate& parent1, const Candidate& parent2) const
    {
        /* Edge case. No point in performing the mutations if the parents are the same. */
        if (parent1 == parent2) return std::make_pair(parent1, parent2);

        switch (crossover_method_)
        {
            case CrossoverMethod::single_point:
                return nPointCrossover(parent1, parent2, crossover_rate_, 1);
            case CrossoverMethod::two_point:
                return nPointCrossover(parent1, parent2, crossover_rate_, 2);
            case CrossoverMethod::n_point:
                return nPointCrossover(parent1, parent2, crossover_rate_, num_crossover_points_);
            case CrossoverMethod::uniform:
                return uniformCrossover(parent1, parent2, crossover_rate_);
            case CrossoverMethod::custom:
                return customCrossover(parent1, parent2, crossover_rate_);
            default:
                assert(false);	/* Invalid crossover method. Shouldn't get here. */
                std::abort();
        }
    }

    inline void IntegerGA::mutate(Candidate& child) const
    {
        switch (mutation_method_)
        {
            case MutationMethod::standard:
                standardMutate(child, mutation_rate_, swap_rate_, inversion_rate_, base_);
                break;
            case MutationMethod::custom:
                customMutate(child, mutation_rate_);
                break;
            default:
                assert(false);	/* Invalid mutation method. Shouldn't get here. */
                std::abort();
        }
    }

    inline IntegerGA::CandidatePair IntegerGA::nPointCrossover(const Candidate& parent1, const Candidate& parent2, double pc, size_t n)
    {
        using namespace std;
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(0.0 <= pc && pc <= 1.0);

        Candidate child1(parent1), child2(parent2);

        /* Perform crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            /* Generate n (or less, but at least 1) number of random loci. */
            unordered_set<size_t> loci;
            for (size_t i = 0; i < n; i++)
            {
                loci.insert(rng::randomInt(size_t{ 1 }, parent1.chromosome.size() - 1));
            }

            /* Count how many loci are after each gene. */
            vector<size_t> loci_after;
            loci_after.reserve(parent1.chromosome.size());

            size_t loci_left = loci.size();
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                if ((loci_left > 0) && loci.contains(i)) loci_left--;
                loci_after.push_back(loci_left);
            }

            /* Perform the crossover. */
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* Swap the bits if there are an odd number of loci after the gene. */
                if (loci_after[i] % 2)
                {
                    child1.chromosome[i] = parent2.chromosome[i];
                    child2.chromosome[i] = parent1.chromosome[i];
                }
            }
            /* Check if the children will need evaluation. */
            if (child1 != parent1)
            {
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }
        }

        return make_pair(child1, child2);
    }

    inline IntegerGA::CandidatePair IntegerGA::uniformCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
    {
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(0.0 <= pc && pc <= 1.0);

        Candidate child1(parent1), child2(parent2);

        /* Perform crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* Swap each bit with 0.5 probability. */
                if (rng::randomBool())
                {
                    child1.chromosome[i] = parent2.chromosome[i];
                    child2.chromosome[i] = parent1.chromosome[i];
                }
            }
            /* Check if the children will need evaluation. */
            if (child1 != parent1)
            {
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }
        }

        return std::make_pair(child1, child2);
    }

    inline void IntegerGA::standardMutate(Candidate& child, double pm, double ps, double pi, size_t base_)
    {
        assert(0.0 <= pm && pm <= 1.0);
        assert(0.0 <= ps && ps <= 1.0);
        assert(0.0 <= pi && pi <= 1.0);
        assert(base_ > 1);

        /* Calc number of mutated genes. */
        double mean = child.chromosome.size() * pm;
        double SD = child.chromosome.size() * pm * (1.0 - pm);

        size_t mutation_count = size_t(std::round(rng::randomNormal(mean, SD)));
        mutation_count = std::clamp(mutation_count, size_t{ 0 }, child.chromosome.size());

        /* The child will (very likely) be changed, and will need to be evaluated. */
        if (mutation_count > 0) child.is_evaluated = false;

        /*
        * Change mutation_count number of random genes.
        * One gene may be changed multiple times, but it's uncommon for long chromosomes, and not important.
        */
        for (size_t i = 0; i < mutation_count; i++)
        {
            size_t idx = rng::randomIdx(child.chromosome.size());
            child.chromosome[idx] = rng::randomInt(size_t{ 0 }, base_ - 1);
        }

        /* Perform swap with ps probability. */
        if (rng::randomReal() <= ps)
        {
            /* r1 and r2 might be the same index, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(child.chromosome.size());
            size_t r2 = rng::randomIdx(child.chromosome.size());
            std::swap(child.chromosome[r1], child.chromosome[r2]);

            if (child.chromosome[r1] != child.chromosome[r2]) child.is_evaluated = false;
        }

        /* Perform inversion with pi probability. */
        if (rng::randomReal() <= pi)
        {
            /* Pick a random range of genes. The bounds of the range may be the same, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(child.chromosome.size());
            size_t r2 = rng::randomIdx(child.chromosome.size());
            auto [idx1, idx2] = std::minmax(r1, r2);

            std::reverse(child.chromosome.begin() + idx1, child.chromosome.begin() + idx2 + 1);

            if (r1 != r2) child.is_evaluated = false;
        }
    }

} // namespace genetic_algorithm

#endif // !GA_INTEGER_GA_H