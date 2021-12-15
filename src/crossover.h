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
* This file contains all the crossover operators for the algorithms.
* Consists of the base Crossover classes that provide the interface for all the crossover operators,
* and implementations of some commonly used crossover operators for all of the encoding types used
* in the algorithms (binary, real, permutation, integer).
*
* @file crossover.h
*/

#ifndef GA_CROSSOVER_H
#define GA_CROSSOVER_H


#include <utility>
#include <cstddef>
#include "candidate.h"

/** Crossover operators for the algorithms. */
namespace genetic_algorithm::crossover
{
    /**
    * Base class used for all of the crossovers. \n
    * Every crossover operator should take 2 Candidates (parents), and return 2 children.
    */
    template<regular_hashable GeneType>
    class Crossover
    {
    public:
        /** Construct a crossover operator with the default parameters. */
        Crossover() = default;

        /**
        * Construct a crossover operator with @p pc as the crossover rate.
        * 
        * @param pc The crossover probability.
        */
        explicit Crossover(double pc);

        virtual ~Crossover() = default;

        /**
        * Sets the crossover rate used in the algorithm to @p pc. \n
        * The crossover rate must be on the closed interval [0.0, 1.0].
        *
        * @param pc The probability of crossover being performed on 2 selected individuals.
        */
        void crossover_rate(double pc);

        /** @returns The crossover rate set for the crossover operator. */
        [[nodiscard]]
        double crossover_rate() const noexcept { return pc_; }

        /**
        * Perform the crossover operator on 2 individuals with the set probability.
        * 
        * @param parent1 The first parent.
        * @param parent2 The second parent.
        * @returns The pair of children resulting from the crossover.
        */
        CandidatePair<GeneType> operator()(const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const;

    protected:

        /* The actual crossover function. Performs the crossover every time and does nothing else. */
        virtual CandidatePair<GeneType> crossover(const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const = 0;

        double pc_ = 0.8;   /* The crossover rate used with the crossover operator. */
    };

    /** 
    * Base class used for the crossover operators where each gene has a lower and upper bound.
    * Used for the real coded crossover operators.
    */
    template<regular_hashable GeneType>
    class BoundedCrossover : public Crossover<GeneType>
    {
    public:
        /**
        * Construct a bounded crossover operator with @p bounds used as the limits for the genes.
        * 
        * @param bounds The bounds of the genes.
        */
        BoundedCrossover(const std::vector<std::pair<GeneType, GeneType>>& bounds);

        /**
        * Construct a bounded crossover operator with @p bounds used as the limits for the genes,
        * and @p pc as the crossover probability.
        */
        BoundedCrossover(double pc, const std::vector<std::pair<GeneType, GeneType>>& bounds);

        virtual ~BoundedCrossover() = default;

        /**
        * Sets the boundaries of the genes. \n
        * Each element of the @p limits vector must contain the lower and upper bounds of the corresponding gene.
        * (The min and max values of the gene.) \n
        * 
        * The number of elements in @p limits must be the same as the length of the chromosomes used,
        * and the lower bound must not be higher than the upper bound for any of the genes. \n
        *
        * @param limits The bounds of the genes.
        */
        void limits(const std::vector<std::pair<GeneType, GeneType>>& limits);

        /** @returns The bounds set for the crossover operator. */
        [[nodiscard]]
        std::vector<std::pair<GeneType, GeneType>> limits() const noexcept { return bounds_; }

    protected:
        std::vector<std::pair<GeneType, GeneType>> bounds_;     /* The lower and upper bounds for each gene of the chromosome. */
    };

    /** Crossover operator templates for any gene type. */
    namespace generic
    {
        /* General n-point crossover template for any gene type. */
        template<regular_hashable GeneType, size_t N>
        class NPointCrossoverTempl : public Crossover<GeneType>
        {
        public:
            using Crossover<GeneType>::Crossover;
        private:
            CandidatePair<GeneType> crossover(const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
        };

        /* General uniform crossover template for any gene type. */
        template<regular_hashable GeneType>
        class UniformCrossoverTempl : public Crossover<GeneType>
        {
        public:
            using Crossover<GeneType>::Crossover;
        private:
            CandidatePair<GeneType> crossover(const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;
        };
    }

    /** Predefined crossover operators for the binary encoded genetic algorithms (BinaryGA). */
    namespace binary
    {
        /**
        * Standard single-point crossover operator for the binary encoded algorithms. \n
        * A random crossover point (locus) is selected and the genes before the locus are swapped
        * between the parents to create the children.
        */
        using SinglePoint = generic::NPointCrossoverTempl<char, 1U>;

        /**
        * Two-point crossover operator for the binary encoded algorithms. \n
        * Two random crossover points are selected, and the genes between the 2 point are
        * swapped between the parents in order to create the children. \n
        * (Works as if 2 consecutive single-point crossovers were performed on the parents.)
        */
        using TwoPoint = generic::NPointCrossoverTempl<char, 2U>;

        /**
        * Uniform crossover operator for the binary encoded algorithms. \n
        * Every gene of the chromosomes is swapped with 0.5 probability between the parents to create the children.
        */
        using Uniform = generic::UniformCrossoverTempl<char>;
    }

    /** Predefined crossover operators for the real encoded genetic algorithms (RCGA). */
    namespace real
    {
        /**
        * Arithmetic crossover operator for the RCGA.
        */
        class Arithmetic : public BoundedCrossover<double>
        {
        public:
            using BoundedCrossover::BoundedCrossover;
        private:
            CandidatePair<double> crossover(const Candidate<double>& parent1, const Candidate<double>& parent2) const override;
        };

        /**
        * BLX-Alpha (blend) crossover operator for the RCGA.
        */
        class BLXa : public BoundedCrossover<double>
        {
        public:
            /**
            * Construct a BLX-alpha crossover operator with @p limits used as the bounds of the genes.
            * 
            * @param limits The lower and upper bounds of each gene.
            */
            BLXa(const std::vector<std::pair<double, double>>& limits);

            /**
            * Construct a BLX-alpha crossover operator with the specified parameters.
            * 
            * @param pc The crossover probability used.
            * @param limits The lower and upper bounds of each gene.
            * @param alpha The alpha parameter of the crossover.
            */
            BLXa(double pc, const std::vector<std::pair<double, double>>& limits, double alpha = 0.5);

            /**
            * Sets the alpha parameter for the crossover. \n
            * This parameter controls the length of the interval the children's genes are randomly chosen from
            * (larger alpha -> larger interval). \n
            * Must be nonnegative, with the ideal value being around 0.5.
            *
            * @param alpha The alpha parameter of the BLX-alpha crossover.
            */
            void alpha(double alpha);

            /** @returns The alpha parameter set for the operator. */
            [[nodiscard]]
            double alpha() const noexcept { return alpha_; }

        private:
            CandidatePair<double> crossover(const Candidate<double>& parent1, const Candidate<double>& parent2) const override;

            double alpha_ = 0.5;        /* The range parameter (alpha) of the BLX-alpha crossover. */
        };

        /**
        * Simulated binary crossover operator for the RCGA.
        */
        class SimulatedBinary : public BoundedCrossover<double>
        {
        public:
            /**
            * Construct a simulated binary crossover operator with @p limits used as the bounds of the genes.
            *
            * @param limits The lower and upper bounds of each gene.
            */
            SimulatedBinary(const std::vector<std::pair<double, double>>& limits);

            /**
            * Construct a simulated binary crossover operator with the specified parameters.
            * 
            * @param pc The crossover probability used.
            * @param limits The lower and upper bounds of each gene.
            * @param eta
            */
            SimulatedBinary(double pc, const std::vector<std::pair<double, double>>& limits, double eta = 4.0);

            /**
            * Sets the shape parameter (eta) of the simulated binary crossover method. \n
            * Must be nonnegative, typical values are 1-5. \n
            * Larger values will lead to the children's genes to be closer to the parents' genes.
            * 
            * @param eta The eta parameter of the crossover.
            */
            void eta(double eta);

            /** @returns The eta parameter set for the operator. */
            [[nodiscard]]
            double eta() const noexcept { return eta_; }

        private:
            CandidatePair<double> crossover(const Candidate<double>& parent1, const Candidate<double>& parent2) const override;

            double eta_ = 4.0;          /* The shape parameter (eta) of the simulated binary crossover. */
        };

        /**
        * Wright's heuristic crossover operator for the RCGA.
        */
        class Wright : public BoundedCrossover<double>
        {
        public:
            using BoundedCrossover::BoundedCrossover;
        private:
            CandidatePair<double> crossover(const Candidate<double>& parent1, const Candidate<double>& parent2) const override;
        };
    }

    /** Predefined crossover operators for the permutation GAs (PermutationGA). */
    namespace perm
    {
        /**
        * Order crossover operator (OX1) for the PermutationGA.
        * Uses no parameters. Fastest method.
        */
        class Order : public Crossover<size_t>
        {
        public:
            using Crossover::Crossover;
        private:
            CandidatePair<size_t> crossover(const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const override;
        };

        /**
        * Cycle crossover operator (CX) for the PermutationGA.
        * Uses no parameters.
        */
        class Cycle : public Crossover<size_t>
        {
        public:
            using Crossover::Crossover;
        private:
            CandidatePair<size_t> crossover(const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const override;
        };

        /**
        * Edge assembly crossover operator (EAX) for the PermutationGA.
        * Uses no parameters. Slowest method.
        */
        class Edge : public Crossover<size_t>
        {
        public:
            using Crossover::Crossover;
        private:
            CandidatePair<size_t> crossover(const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const override;
        };

        /**
        * Partially mapped crossover operator (PMX) for the PermutationGA.
        * Uses no parameters.
        */
        class PMX : public Crossover<size_t>
        {
        public:
            using Crossover::Crossover;
        private:
            CandidatePair<size_t> crossover(const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const override;
        };
    }

    /** Predefined crossover operators for the integer encoded genetic algorithms (IntegerGA). */
    namespace integer
    {
        /**
        * Standard single-point crossover operator for the integer encoded algorithms. \n
        * A random crossover point (locus) is selected and the genes before the locus are swapped
        * between the parents to create the children.
        */
        using SinglePoint = generic::NPointCrossoverTempl<size_t, 1U>;

        /**
        * Two-point crossover operator for the integer encoded algorithms. \n
        * Two random crossover points are selected, and the genes between the 2 point are
        * swapped between the parents in order to create the children. \n
        * (Works as if 2 consecutive single-point crossovers were performed on the parents.)
        */
        using TwoPoint = generic::NPointCrossoverTempl<size_t, 2U>;

        /**
        * Uniform crossover operator for the integer encoded algorithms. \n
        * Every gene of the chromosomes is swapped with 0.5 probability between the parents to create the children.
        */
        using Uniform = generic::UniformCrossoverTempl<size_t>;
    }

} // namespace genetic_algorithm::crossover


/* IMPLEMENTATION */

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <limits>
#include <stdexcept>
#include <cassert>
#include "rng.h"
#include "mo_detail.h"

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
    inline CandidatePair<T> Crossover<T>::operator()(const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        assert(0.0 <= pc_ && pc_ <= 1.0);

        /* Only need to perform the crossover with the set pc probability. Return early with (1 - pc) probability. */
        if (rng::randomReal() >= pc_)
        {
            return { parent1, parent2 };
        }

        /* 
        * If the parents are the same, the crossover doesn't need to be performed.
        * This assumes that with 2 parents that have the same chromosomes, the children would be the same
        * as the parents, which is true for every crossover operator implemented, but
        * could be an issue for user defined crossovers.
        */
        if (parent1 == parent2)
        {
            return { parent1, parent2 };
        }

        /* Perform the actual crossover. */
        auto [child1, child2] = crossover(parent1, parent2);

        child1.is_evaluated = false;
        child2.is_evaluated = false;

        /* 
        * Check if either of the children are the same as one of the parents.
        * (This can happen in edge cases even if the parents are different)
        * If one of the children are the same as one of the parents, then the fitness function
        * evaluation for that child can be skipped (if the fitness function is the same) by assigning
        * it the same fitness as the parent.
        */
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

        return { child1, child2 };
    }

    template<regular_hashable GeneType>
    inline BoundedCrossover<GeneType>::BoundedCrossover(const std::vector<std::pair<GeneType, GeneType>>& bounds) 
        : Crossover<GeneType>()
    {
        limits(bounds);
    }

    template<regular_hashable GeneType>
    inline BoundedCrossover<GeneType>::BoundedCrossover(double pc, const std::vector<std::pair<GeneType, GeneType>>& bounds) 
        : Crossover<GeneType>(pc)
    {
        limits(bounds);
    }

    template<regular_hashable GeneType>
    inline void BoundedCrossover<GeneType>::limits(const std::vector<std::pair<GeneType, GeneType>>& limits)
    {
        if (std::any_of(limits.begin(), limits.end(), [](std::pair<double, double> b) {return b.first > b.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

        bounds_ = limits;
    }

    namespace generic
    {
        template<regular_hashable T, size_t N>
        inline CandidatePair<T> NPointCrossoverTempl<T, N>::crossover(const Candidate<T>& parent1, const Candidate<T>& parent2) const
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());

            Candidate<T> child1{ parent1 }, child2{ parent2 };

            /* No reason to have more loci than genes. */
            size_t num_loci = std::min(N, parent1.chromosome.size());
            /* Generate n (or less, but at least 1) number of unique random loci. */
            std::unordered_set<size_t> loci;
            for (size_t i = 0; i < num_loci; i++)
            {
                loci.insert(rng::randomInt(size_t{ 1 }, parent1.chromosome.size() - 1));
            }

            /* Count how many loci are after each gene. */
            std::vector<size_t> loci_after;
            loci_after.reserve(parent1.chromosome.size());

            for (size_t i = 0, loci_left = loci.size(); i < parent1.chromosome.size(); i++)
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

        template<regular_hashable T>
        inline CandidatePair<T> UniformCrossoverTempl<T>::crossover(const Candidate<T>& parent1, const Candidate<T>& parent2) const
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());

            Candidate child1{ parent1 }, child2{ parent2 };

            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* Swap each gene with 0.5 probability. */
                if (rng::randomBool())
                {
                    child1.chromosome[i] = parent2.chromosome[i];
                    child2.chromosome[i] = parent1.chromosome[i];
                }
            }

            return { child1, child2 };
        }
    }

    namespace real
    {
        inline BLXa::BLXa(const std::vector<std::pair<double, double>>& limits) :
            BoundedCrossover(limits)
        {}

        inline BLXa::BLXa(double pc, const std::vector<std::pair<double, double>>& limits, double alpha) :
            BoundedCrossover(pc, limits)
        {
            this->alpha(alpha);
        }

        inline void BLXa::alpha(double alpha)
        {
            if (!(0.0 <= alpha && alpha <= std::numeric_limits<double>::max()))
            {
                throw std::invalid_argument("Alpha must be a nonnegative, finite value.");
            }

            alpha_ = alpha;
        }

        inline SimulatedBinary::SimulatedBinary(const std::vector<std::pair<double, double>>& limits) :
            BoundedCrossover(limits)
        {}

        inline SimulatedBinary::SimulatedBinary(double pc, const std::vector<std::pair<double, double>>& limits, double eta) :
            BoundedCrossover(pc, limits)
        {
            this->eta(eta);
        }

        inline void SimulatedBinary::eta(double eta)
        {
            if (!(0.0 <= eta && eta <= std::numeric_limits<double>::max()))
            {
                throw std::invalid_argument("Eta must be a nonnegative, finite value.");
            }

            eta_ = eta;
        }

        inline CandidatePair<double> Arithmetic::crossover(const Candidate<double>& parent1, const Candidate<double>& parent2) const
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());

            Candidate child1{ parent1 }, child2{ parent2 };

            double alpha = rng::randomReal();
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                child1.chromosome[i] = alpha * parent1.chromosome[i] + (1.0 - alpha) * parent2.chromosome[i];
                child2.chromosome[i] = (1.0 - alpha) * parent1.chromosome[i] + alpha * parent2.chromosome[i];
            }
            /* No bounds check, the generated children's genes will always be within the bounds if the parents' genes were also within them. */

            return { child1, child2 };
        }

        inline CandidatePair<double> BLXa::crossover(const Candidate<double>& parent1, const Candidate<double>& parent2) const
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(parent1.chromosome.size() == bounds_.size());
            assert(alpha_ >= 0.0);

            Candidate child1{ parent1 }, child2{ parent2 };

            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* Calc interval to generate the childrens genes on. */
                auto [range_min, range_max] = std::minmax(parent1.chromosome[i], parent2.chromosome[i]);
                double range_ext = alpha_ * (range_max - range_min);
                /* Generate genes from an uniform distribution on the interval. */
                child1.chromosome[i] = rng::randomReal(range_min - range_ext, range_max + range_ext);
                child2.chromosome[i] = rng::randomReal(range_min - range_ext, range_max + range_ext);
                /* The children's genes might be outside the allowed interval. */
                child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds_[i].first, bounds_[i].second);
                child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds_[i].first, bounds_[i].second);
            }

            return { child1, child2 };
        }

        inline CandidatePair<double> SimulatedBinary::crossover(const Candidate<double>& parent1, const Candidate<double>& parent2) const
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(parent1.chromosome.size() == bounds_.size());
            assert(eta_ > 0.0);

            Candidate child1{ parent1 }, child2{ parent2 };

            /* Generate beta from the appropriate distribution. */
            double u = rng::randomReal();
            double beta = (u <= 0.5) ? std::pow(2.0 * u, 1.0 / (eta_ + 1.0)) : 
                                       std::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta_ + 1.0));

            /* Perform crossover. */
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                child1.chromosome[i] = 0.5 * ((1 - beta) * parent1.chromosome[i] + (1 + beta) * parent2.chromosome[i]);
                child2.chromosome[i] = 0.5 * ((1 + beta) * parent1.chromosome[i] + (1 - beta) * parent2.chromosome[i]);
                /* The children's genes might be outside the allowed interval. */
                child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds_[i].first, bounds_[i].second);
                child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds_[i].first, bounds_[i].second);
            }

            return { child1, child2 };
        }

        inline CandidatePair<double> Wright::crossover(const Candidate<double>& parent1, const Candidate<double>& parent2) const
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(parent1.chromosome.size() == bounds_.size());

            Candidate child1(parent1), child2(parent2);

            /* p1 is always the better parent. */
            const Candidate<double>* p1 = detail::paretoCompare(parent1.fitness, parent2.fitness) ? &parent2 : &parent1;
            const Candidate<double>* p2 = detail::paretoCompare(parent1.fitness, parent2.fitness) ? &parent1 : &parent2;
            /* Get random weights. */
            double w1 = rng::randomReal();
            double w2 = rng::randomReal();
            /* Perform crossover. */
            for (size_t i = 0; i < p1->chromosome.size(); i++)
            {
                child1.chromosome[i] = w1 * (p1->chromosome[i] - p2->chromosome[i]) + p1->chromosome[i];
                child2.chromosome[i] = w2 * (p1->chromosome[i] - p2->chromosome[i]) + p1->chromosome[i];
                /* The children's genes might be outside the allowed intervals. */
                child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds_[i].first, bounds_[i].second);
                child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds_[i].first, bounds_[i].second);
            }

            return { child1, child2 };
        }

    } // namespace crossover:real

    namespace perm
    {
        inline CandidatePair<size_t> Order::crossover(const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
        {
            using namespace std;
            assert(parent1.chromosome.size() == parent2.chromosome.size());

            /* Pick a random range of genes. */
            size_t r1 = rng::randomIdx(parent1.chromosome.size());
            size_t r2 = rng::randomIdx(parent1.chromosome.size());
            const auto [idx1, idx2] = minmax(r1, r2);

            /* Edge case. The entire chromosomes are swapped. */
            if (idx1 == 0 && idx2 == parent1.chromosome.size() - 1)
            {
                return { parent2, parent1 };
            }

            /* Gather the genes that will go directly from parent1 -> child1, and parent2 -> child2. (Not using the constructor is intentional.) */
            unordered_set<size_t> direct1, direct2;
            for (auto gene = parent1.chromosome.begin() + idx1; gene != parent1.chromosome.begin() + idx2 + 1; gene++) direct1.insert(*gene);
            for (auto gene = parent2.chromosome.begin() + idx1; gene != parent2.chromosome.begin() + idx2 + 1; gene++) direct2.insert(*gene);

            /* Gather the remaining genes (not in the range) from the other parent. */
            vector<size_t> cross1;    /* Segment gathered from parent2 -> child1. */
            vector<size_t> cross2;    /* Segment gathered from parent1 -> child2. */
            cross1.reserve(parent2.chromosome.size() - idx2 + idx1 - 1);
            cross2.reserve(parent1.chromosome.size() - idx2 + idx1 - 1);
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* If a gene is not taken directly from the corresponding, then take it from the other parent. */
                if (!direct1.contains(parent2.chromosome[i])) cross1.push_back(parent2.chromosome[i]);
                if (!direct2.contains(parent1.chromosome[i])) cross2.push_back(parent1.chromosome[i]);
            }

            /* Construct the children: child1 = direct1 + cross1, child2 = direct2 + cross2 */
            Candidate<size_t> child1, child2;
            child1.chromosome.reserve(parent1.chromosome.size());
            child2.chromosome.reserve(parent1.chromosome.size());
            /* child1 */
            child1.chromosome.insert(child1.chromosome.end(), cross1.begin(), cross1.begin() + idx1);
            child1.chromosome.insert(child1.chromosome.end(), parent1.chromosome.begin() + idx1, parent1.chromosome.begin() + idx2 + 1);
            child1.chromosome.insert(child1.chromosome.end(), cross1.begin() + idx1, cross1.end());
            /* child2 */
            child2.chromosome.insert(child2.chromosome.end(), cross2.begin(), cross2.begin() + idx1);
            child2.chromosome.insert(child2.chromosome.end(), parent2.chromosome.begin() + idx1, parent2.chromosome.begin() + idx2 + 1);
            child2.chromosome.insert(child2.chromosome.end(), cross2.begin() + idx1, cross2.end());

            return { child1, child2 };
        }

        inline CandidatePair<size_t> Cycle::crossover(const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
        {
            using namespace std;
            assert(parent1.chromosome.size() == parent2.chromosome.size());

            /* Identify all cycles. */
            vector<vector<size_t>> cycles;
            /* Copies of parent chromosomes so they can be changed without issues. */
            vector<size_t> chrom1 = parent1.chromosome;
            vector<size_t> chrom2 = parent2.chromosome;
            while (!chrom1.empty())
            {
                vector<size_t> cycle;
                /* Always start the cycle at chrom1[0]. Old cycles are removed the chromosomes. */
                size_t cur_pos = 0;
                size_t top_val = chrom1[cur_pos];
                cycle.push_back(top_val);

                while (chrom2[cur_pos] != chrom1[0])
                {
                    /* Look for the bottom value at this pos (chrom2[pos]) in chrom1. This is the new pos and top value. */
                    cur_pos = static_cast<size_t>(find(chrom1.begin(), chrom1.end(), chrom2[cur_pos]) - chrom1.begin());
                    top_val = chrom1[cur_pos];
                    /* Add the new top value to the cycle. */
                    cycle.push_back(top_val);
                    /* Keep going until the bottom value at this pos isn't the cycle start value. (cycle complete) */
                }

                /* Delete the values in this cycle from chrom1 and chrom2. (Without changing the order of the remaining genes.) */
                for (size_t i = 0; i < cycle.size(); i++)
                {
                    chrom1.erase(find(chrom1.begin(), chrom1.end(), cycle[i]));
                    chrom2.erase(find(chrom2.begin(), chrom2.end(), cycle[i]));
                }
                /* Add this cycle to the cycles. */
                cycles.push_back(cycle);
            }

            /* Construct the children from the cycles. */
            Candidate child1{ parent1 }, child2{ parent2 };

            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* Find which cycle has the gene. */
                size_t cycle_num = 0;
                for (size_t j = 0; j < cycles.size(); j++)
                {
                    if (find(cycles[j].begin(), cycles[j].end(), parent1.chromosome[i]) != cycles[j].end())
                    {
                        cycle_num = j + 1;
                        break;
                    }
                }
                /* Even cycle genes are swapped parent1->child2 and parent2->child1. */
                if (cycle_num % 2 == 0)
                {
                    child1.chromosome[i] = parent2.chromosome[i];
                    child2.chromosome[i] = parent1.chromosome[i];
                }
                /* Odd cycle genes were already handled when initializing the children. */
            }

            return { child1, child2 };
        }

        inline CandidatePair<size_t> Edge::crossover(const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
        {
            using namespace std;
            using NList = vector<unordered_set<size_t>>;

            assert(parent1.chromosome.size() == parent2.chromosome.size());

            Candidate<size_t> child1, child2;
            child1.chromosome.reserve(parent1.chromosome.size());
            child2.chromosome.reserve(parent2.chromosome.size());

            /* Construct neighbour list based on parents. The first and last genes are not neighbours. */
            size_t len = parent1.chromosome.size();
            NList nl1(len);
            /* Neighbours of first and last genes. */
            nl1[parent1.chromosome.front()] = { parent1.chromosome[1] };
            nl1[parent1.chromosome.back()] = { parent1.chromosome[len - 2] };
            nl1[parent2.chromosome.front()].insert(parent2.chromosome[1]);
            nl1[parent2.chromosome.back()].insert(parent2.chromosome[len - 2]);
            /* Neighbours of all other genes. */
            for (size_t i = 1; i < len - 1; i++)
            {
                nl1[parent1.chromosome[i]].insert(parent1.chromosome[i + 1]);
                nl1[parent1.chromosome[i]].insert(parent1.chromosome[i - 1]);

                nl1[parent2.chromosome[i]].insert(parent2.chromosome[i + 1]);
                nl1[parent2.chromosome[i]].insert(parent2.chromosome[i - 1]);
            }
            NList nl2 = nl1;    /* Copy for child2. */

            /* Generate child1. */
            size_t gene = parent1.chromosome[0];
            vector<size_t> remaining_genes = parent1.chromosome;

            while (child1.chromosome.size() != parent1.chromosome.size())
            {
                /* Append gene to the child, and remove it from all neighbour lists. */
                child1.chromosome.push_back(gene);
                remaining_genes.erase(remove(remaining_genes.begin(), remaining_genes.end(), gene), remaining_genes.end());
                for (auto& neighbours : nl1)
                {
                    neighbours.erase(gene);
                }
                if (child1.chromosome.size() == parent1.chromosome.size()) break;

                /* Determine next gene that will be added to the child. */
                /* If gene's neighbour list is empty, gene = random node not already in child. */
                if (nl1[gene].empty())
                {
                    gene = remaining_genes[rng::randomIdx(remaining_genes.size())];
                }
                else /* gene's neighbour list is not empty, gene = neighbour of gene with fewest neighbours (random if tie). */
                {
                    /* Find gene's neighbour with fewest neighbours. */
                    size_t nb = *min_element(nl1[gene].begin(), nl1[gene].end(),
                    [&nl1](const size_t& lhs, const size_t& rhs)
                    {
                        return (nl1[lhs].size() < nl1[rhs].size());
                    });
                    size_t min_neighbour_count = nl1[nb].size();

                    /* Determine possible nodes (neighbours of gene with min_neighbour_count neighbours). */
                    vector<size_t> possible_nodes;
                    for (const auto& neighbour : nl1[gene])
                    {
                        if (nl1[neighbour].size() == min_neighbour_count)
                        {
                            possible_nodes.push_back(neighbour);
                        }
                    }

                    gene = possible_nodes[rng::randomIdx(possible_nodes.size())];
                }
            }

            /* Same process to get child2. */
            gene = parent2.chromosome[0];
            remaining_genes = parent2.chromosome;
            while (child2.chromosome.size() != parent2.chromosome.size())
            {
                child2.chromosome.push_back(gene);
                remaining_genes.erase(remove(remaining_genes.begin(), remaining_genes.end(), gene), remaining_genes.end());
                for (auto& neighbours : nl2)
                {
                    neighbours.erase(gene);
                }
                if (child2.chromosome.size() == parent2.chromosome.size()) break;

                if (nl2[gene].empty())
                {
                    gene = remaining_genes[rng::randomIdx(remaining_genes.size())];
                }
                else
                {
                    size_t nb = *min_element(nl2[gene].begin(), nl2[gene].end(),
                    [&nl2](const size_t& lhs, const size_t& rhs)
                    {
                        return (nl2[lhs].size() < nl2[rhs].size());
                    });
                    size_t min_neighbour_count = nl2[nb].size();

                    vector<size_t> possible_nodes;
                    for (const auto& neighbour : nl2[gene])
                    {
                        if (nl2[neighbour].size() == min_neighbour_count)
                        {
                            possible_nodes.push_back(neighbour);
                        }
                    }

                    gene = possible_nodes[rng::randomIdx(possible_nodes.size())];
                }
            }

            return { child1, child2 };
        }

        inline CandidatePair<size_t> PMX::crossover(const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const
        {
            using namespace std;
            assert(parent1.chromosome.size() == parent2.chromosome.size());

            /* Init children so the last step of the crossover can be skipped. */
            Candidate child1{ parent2 }, child2{ parent1 };

            /* Pick a random range of genes. The bounds of the range may be the same, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(parent1.chromosome.size());
            size_t r2 = rng::randomIdx(parent1.chromosome.size());
            const auto [idx1, idx2] = minmax(r1, r2);

            /* Edge case. The entire chromosomes are copied directly. */
            if (idx1 == 0 && idx2 == parent1.chromosome.size() - 1)
            {
                return { parent1, parent2 };
            }

            /* Copy genes in the range from the corresponding parent. */
            for (size_t i = idx1; i <= idx2; i++)
            {
                child1.chromosome[i] = parent1.chromosome[i];
                child2.chromosome[i] = parent2.chromosome[i];
            }
            /* Note ranges that were copied from parents to check if they contain an element. (Not using the constructor is intentional.) */
            unordered_set<size_t> p1_range, p2_range;
            for (auto gene = parent1.chromosome.begin() + idx1; gene != parent1.chromosome.begin() + idx2 + 1; gene++) p1_range.insert(*gene);
            for (auto gene = parent2.chromosome.begin() + idx1; gene != parent2.chromosome.begin() + idx2 + 1; gene++) p2_range.insert(*gene);

            /* Get the rest of the children's genes from the other parents. */
            for (size_t i = idx1; i <= idx2; i++)
            {
                /* child1 */
                if (!p1_range.contains(parent2.chromosome[i]))
                {
                    size_t cur_pos = i;
                    while (idx1 <= cur_pos && cur_pos <= idx2)
                    {
                        size_t val = parent1.chromosome[cur_pos];
                        cur_pos = static_cast<size_t>(find(parent2.chromosome.begin(), parent2.chromosome.end(), val) - parent2.chromosome.begin());
                    }
                    child1.chromosome[cur_pos] = parent2.chromosome[i];
                }
                /* child2 */
                if (!p2_range.contains(parent1.chromosome[i]))
                {
                    size_t cur_pos = i;
                    while (idx1 <= cur_pos && cur_pos <= idx2)
                    {
                        size_t val = parent2.chromosome[cur_pos];
                        cur_pos = static_cast<size_t>(find(parent1.chromosome.begin(), parent1.chromosome.end(), val) - parent1.chromosome.begin());
                    }
                    child2.chromosome[cur_pos] = parent1.chromosome[i];
                }
            }
            /* Copy any not yet in child positions to the children from the other parents. (Already done at the initialization of the children.) */

            return { child1, child2 };
        }

    } // namespace crossover::perm

} // namespace crossover

#endif // !GA_CROSSOVER_H