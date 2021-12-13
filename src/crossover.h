/*
*  MIT License
*
*  Copyright (c) 2021 Kriszti�n Rug�si
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
* This file contains all the crossover operators for the algorithms.
* Consists of the base Crossover class that provides the interface for all the crossover operators,
* and implementations of commonly used crossover operators for binary, real, permutation, integer genetic
* algorithms.
*
* @file crossover.h
*/

#ifndef GA_CROSSOVER_H
#define GA_CROSSOVER_H


#include <utility>
#include <cstddef>
#include "candidate.h"

namespace genetic_algorithm::crossover
{
    template<regular_hashable GeneType>
    struct Crossover
    {
        Crossover() = default;

        explicit Crossover(double pc);

        void crossover_rate(double pc);

        [[nodiscard]] double crossover_rate() const noexcept { return pc_; }

        std::pair<Candidate<GeneType>, Candidate<GeneType>> operator()(const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2);

    protected:

        virtual std::pair<Candidate<GeneType>, Candidate<GeneType>> crossover(const Candidate<GeneType>& parent1,
                                                                              const Candidate<GeneType>& parent2) const = 0;

        double pc_ = 0.8;   /* The crossover rate used with the crossover operator. */
    };

    namespace binary
    {
        struct SinglePoint final : Crossover<char>
        {
            using Crossover::Crossover;

        private:           
            std::pair<Candidate<char>, Candidate<char>> crossover(const Candidate<char>& parent1, const Candidate<char>& parent2) const override;
        };

        struct TwoPoint final : Crossover<char>
        {
            using Crossover::Crossover;

        private:
            std::pair<Candidate<char>, Candidate<char>> crossover(const Candidate<char>& parent1, const Candidate<char>& parent2) const override;
        };

        struct NPoint final : Crossover<char>
        {
            NPoint() = default;
            NPoint(double pc, size_t n);

            void num_crossover_points(size_t n);
            [[nodiscard]] size_t num_crossover_points() const noexcept { return n_; }


        private:
            size_t n_ = 1;  /* The number of crossover points used. */

            std::pair<Candidate<char>, Candidate<char>> crossover(const Candidate<char>& parent1, const Candidate<char>& parent2) const override;
        };

        struct Uniform final : Crossover<char>
        {
            using Crossover::Crossover;

        private:
            std::pair<Candidate<char>, Candidate<char>> crossover(const Candidate<char>& parent1, const Candidate<char>& parent2) const override;
        };
    }
    namespace real
    {

    }
    namespace perm
    {

    }
    namespace integer
    {

    }

} // namespace genetic_algorithm::crossover

namespace genetic_algorithm::detail
{
    template<typename T>
    inline std::pair<Candidate<T>, Candidate<T>> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t n);
}


/* IMPLEMENTATION */

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <cassert>
#include "rng.h"

namespace genetic_algorithm::crossover
{
    template<regular_hashable GeneType>
    inline Crossover<GeneType>::Crossover(double pc)
    {
        crossover_rate(pc);
    }

    template<regular_hashable GeneType>
    inline void Crossover<GeneType>::crossover_rate(double pc)
    {
        if (!(0.0 <= pc && pc <= 1.0))
        {
            throw std::invalid_argument("The crossover probability must be in the closed range [0.0, 1.0].");
        }

        pc_ = pc;
    }

    template<regular_hashable T>
    inline std::pair<Candidate<T>, Candidate<T>> Crossover<T>::operator()(const Candidate<T>& parent1, const Candidate<T>& parent2)
    {
        /* If the parents are the same, the crossover doesn't need to be performed. */
        if (parent1 == parent2)
        {
            return { parent1, parent2 };
        }
        /* Only need to perform the crossover with the set pc probability. (Return early with (1 - pc) probability.) */
        if (rng::randomReal() >= pc_)
        {
            return { parent1, parent2 };
        }

        /* Perform the actual crossover. */
        auto [child1, child2] = crossover(parent1, parent2);

        child1.is_evaluated = false;
        child2.is_evaluated = false;

        /* Check if either of the children are the same as one of the parents, and won't need to be evaluated. */
        if (child1 == parent1)
        {
            child1.fitness = parent1.fitness;
            child1.is_evaluated = true;
        }
        else if (child1 == parent2)
        {
            child1.fitness = parent2.fitness;
            child1.is_evaluated = true;
        }
        if (child2 == parent1)
        {
            child2.fitness = parent1.fitness;
            child2.is_evaluated = true;
        }
        else if (child2 == parent2)
        {
            child2.fitness = parent2.fitness;
            child2.is_evaluated = true;
        }

        return std::make_pair(child1, child2);
    }

    namespace binary
    {
        inline std::pair<Candidate<char>, Candidate<char>> SinglePoint::crossover(const Candidate<char>& parent1, const Candidate<char>& parent2) const
        {
            return detail::nPointCrossoverImpl(parent1, parent2, 1);
        }

        inline std::pair<Candidate<char>, Candidate<char>> TwoPoint::crossover(const Candidate<char>& parent1, const Candidate<char>& parent2) const
        {
            return detail::nPointCrossoverImpl(parent1, parent2, 2);
        }

        inline std::pair<Candidate<char>, Candidate<char>> NPoint::crossover(const Candidate<char>& parent1, const Candidate<char>& parent2) const
        {
            return detail::nPointCrossoverImpl(parent1, parent2, n_);
        }

        inline NPoint::NPoint(double pc, size_t n)
        {
            crossover_rate(pc);
            num_crossover_points(n);
        }

        inline void NPoint::num_crossover_points(size_t n)
        {
            if (n == 0)
            {
                throw std::invalid_argument("The number of crossover points must be at least 1.");
            }

            n_ = n;
        }

        inline std::pair<Candidate<char>, Candidate<char>> Uniform::crossover(const Candidate<char>& parent1, const Candidate<char>& parent2) const
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());

            Candidate<char> child1(parent1), child2(parent2);

            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* Swap each gene with 0.5 probability. */
                if (rng::randomBool())
                {
                    child1.chromosome[i] = parent2.chromosome[i];
                    child2.chromosome[i] = parent1.chromosome[i];
                }
            }

            return std::make_pair(child1, child2);
        }
    }

    namespace real
    {

    }

    namespace perm
    {

    }

    namespace integer
    {

    }
}

namespace genetic_algorithm::detail
{
    template<typename T>
    inline std::pair<Candidate<T>, Candidate<T>> nPointCrossoverImpl(const Candidate<T>& parent1, const Candidate<T>& parent2, size_t n)
    {
        assert(parent1.chromosome.size() == parent2.chromosome.size());

        Candidate<T> child1(parent1), child2(parent2);

        /* Generate n (or less, but at least 1) number of unique random loci. */
        std::unordered_set<size_t> loci;
        for (size_t i = 0; i < std::min(n, parent1.chromosome.size()); i++)
        {
            loci.insert(rng::randomInt(size_t{ 1 }, parent1.chromosome.size() - 1));
        }

        /* Count how many loci are after each gene. */
        std::vector<size_t> loci_after;
        loci_after.reserve(parent1.chromosome.size());

        size_t loci_left = loci.size();
        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            if (loci_left > 0 && loci.contains(i)) loci_left--;
            loci_after.push_back(loci_left);
        }

        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            /* Swap the genes if there are an odd number of loci after the gene. */
            if (loci_after[i] % 2)
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
        }

        return { child1, child2 };
    }

} // namespace genetic_algorithm::detail

#endif // !GA_CROSSOVER_H