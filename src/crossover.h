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
    namespace binary
    {
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
} // namespace genetic_algorithm::detail

#endif // !GA_CROSSOVER_H