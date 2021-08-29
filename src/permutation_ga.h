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
* This file contains the permutations GA class.
*
* @file permutation_ga.h
*/

#ifndef GA_PERMUTATION_GA_H
#define GA_PERMUTATION_GA_H

#include <cstddef>

#include "base_ga.h"

namespace genetic_algorithm
{
    /**
    * Genetic algorithm that uses permutational encoding. \n
    * The genes of the chromosomes are all unique unsigned integers on [0, chrom_len-1].
    */
    class PermutationGA : public GA<size_t>
    {
    public:

        /**
        * Possible crossover methods that can be used in the PermutationGA. \n
        * Includes some commonly used crossover operators in permutations GAs, but a custom function can also be used
        * to perform the crossovers with the custom option. \n
        * Set the crossover method used in the algorithm with @ref crossover_method.
        */
        enum class CrossoverMethod
        {
            order,	/**< Order crossover operator (OX1). Uses no parameters. Fastest method. */
            cycle,	/**< Cycle crossover operator (CX). Uses no parameters. */
            edge,	/**< Edge assembly crossover operator (EAX). Uses no parameters. Slowest method. */
            pmx,	/**< Partially mapped crossover operator (PMX). Uses no parameters. */
            custom	/**< Custom crossover function defined by the user. @see setCrossoverFunction */
        };

        /**
        * Possible mutation methods that can be used in the PermutationGA. \n
        * Includes commonly used mutation operators in permutations GAs, but a custom mutation function can
        * also be used to perform the mutations with the custom option. \n
        * Set the mutation method used in the algorithm with @ref mutation_method.
        */
        enum class MutationMethod
        {
            swap,		/**< Single-swap mutation operator. Uses no parameters. */
            scramble,	/**< Scramble mutation operator. Uses no parameters. */
            inversion,	/**< Inversion mutation operator. Uses no parameters. */
            custom		/**< Custom mutation function defined by the user. @see setMutationFunction */
        };

        /**
        * Basic contructor for the PermutationGA.
        *
        * @param chrom_len The number of genes in the chromosomes.
        * @param fitness_function The fitness function used in the algorithm to find the maximum of.
        */
        PermutationGA(size_t chrom_len, fitnessFunction_t fitnessFunction);

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

    private:

        CrossoverMethod crossover_method_ = CrossoverMethod::order;
        MutationMethod mutation_method_ = MutationMethod::inversion;

        Candidate generateCandidate() const override;
        CandidatePair crossover(const Candidate& parent1, const Candidate& parent2) const override;
        void mutate(Candidate& child) const override;

        static CandidatePair orderCrossover(const Candidate& parent1, const Candidate& parent2, double pc);
        static CandidatePair cycleCrossover(const Candidate& parent1, const Candidate& parent2, double pc);
        static CandidatePair edgeCrossover(const Candidate& parent1, const Candidate& parent2, double pc);
        static CandidatePair pmxCrossover(const Candidate& parent1, const Candidate& parent2, double pc);

        static void swapMutate(Candidate& child, double pm);
        static void scrambleMutate(Candidate& child, double pm);
        static void inversionMutate(Candidate& child, double pm);
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <algorithm>
#include <random>
#include <vector>
#include <unordered_set>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <cassert>
#include <cstdlib>

#include "rng.h"

namespace genetic_algorithm
{
    inline PermutationGA::PermutationGA(size_t chrom_len, fitnessFunction_t fitnessFunction)
        : GA(chrom_len, fitnessFunction)
    {
    }


    inline void PermutationGA::crossover_method(crossoverFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        crossover_method_ = CrossoverMethod::custom;
        customCrossover = f;
    }

    inline void PermutationGA::crossover_method(CrossoverMethod method)
    {
        if (static_cast<size_t>(method) > 4) throw std::invalid_argument("Invalid crossover method selected.");

        crossover_method_ = method;
    }

    inline PermutationGA::CrossoverMethod PermutationGA::crossover_method() const
    {
        return crossover_method_;
    }


    inline void PermutationGA::mutation_method(mutationFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        mutation_method_ = MutationMethod::custom;
        customMutate = f;
    }


    inline void PermutationGA::mutation_method(MutationMethod method)
    {
        if (static_cast<size_t>(method) > 3) throw std::invalid_argument("Invalid mutation method selected.");

        mutation_method_ = method;
    }

    inline PermutationGA::MutationMethod PermutationGA::mutation_method() const
    {
        return mutation_method_;
    }


    inline PermutationGA::Candidate PermutationGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;

        std::vector<size_t> chrom(chrom_len_);
        std::iota(chrom.begin(), chrom.end(), 0U);
        std::shuffle(chrom.begin(), chrom.end(), rng::prng);

        sol.chromosome = chrom;

        return sol;
    }

    inline PermutationGA::CandidatePair PermutationGA::crossover(const Candidate& parent1, const Candidate& parent2) const
    {
        using namespace std;

        /* Edge case. No point in performing the crossover if the parents are the same. */
        if (parent1 == parent2)
        {
            return make_pair(parent1, parent2);
        }

        Candidate child1, child2;
        switch (crossover_method_)
        {
            case CrossoverMethod::order:
                tie(child1, child2) = orderCrossover(parent1, parent2, crossover_rate_);
                break;
            case CrossoverMethod::cycle:
                tie(child1, child2) = cycleCrossover(parent1, parent2, crossover_rate_);
                break;
            case CrossoverMethod::pmx:
                tie(child1, child2) = pmxCrossover(parent1, parent2, crossover_rate_);
                break;
            case CrossoverMethod::edge:
                tie(child1, child2) = edgeCrossover(parent1, parent2, crossover_rate_);
                break;
            case CrossoverMethod::custom:
                tie(child1, child2) = customCrossover(parent1, parent2, crossover_rate_);
                break;
            default:
                assert(false);	/* Invalid crossover method. Shouldn't get here. */
                abort();
        }

        /* Check if the evaluation of the children can be skipped. */
        /*
        * These checks decrease fitness evals by a lot for short chromosomes:
        *	TSP13 (200pop, 1000gen, 0.9pc):	~168'000 -> ~28'000 fitness evals -(second and last checks added)-> ~20'000 evals
        * Smaller decrease for long chromosomes:
        *	chrom_len=10'000 (50pop, 10gen, 1.0pc): 500 -> ~480-495 fitness evals -(second and last checks added)-> ~475-490 evals
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
        if (child2 == parent2)
        {
            child2.fitness = parent2.fitness;
            child2.is_evaluated = true;
        }
        else if (child2 == parent1)
        {
            child2.fitness = parent1.fitness;
            child2.is_evaluated = true;
        }

        return make_pair(child1, child2);
    }

    inline void PermutationGA::mutate(Candidate& child) const
    {
        switch (mutation_method_)
        {
            case MutationMethod::swap:
                swapMutate(child, mutation_rate_);
                break;
            case MutationMethod::scramble:
                scrambleMutate(child, mutation_rate_);
                break;
            case MutationMethod::inversion:
                inversionMutate(child, mutation_rate_);
                break;
            case MutationMethod::custom:
                customMutate(child, mutation_rate_);
                break;
            default:
                assert(false);	/* Invalid mutation method. Shouldn't get here. */
                std::abort();
        }
    }

    inline PermutationGA::CandidatePair PermutationGA::orderCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
    {
        using namespace std;
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(0.0 <= pc && pc <= 1.0);

        Candidate child1, child2;
        child1.chromosome.reserve(parent1.chromosome.size());
        child2.chromosome.reserve(parent1.chromosome.size());

        /* Perform crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            /* Pick a random range of genes. */
            size_t r1 = rng::randomIdx(parent1.chromosome.size());
            size_t r2 = rng::randomIdx(parent1.chromosome.size());
            const auto [idx1, idx2] = minmax(r1, r2);

            /* Edge case. The entire chromosomes are swapped. */
            if (idx1 == 0 && idx2 == parent1.chromosome.size() - 1) return make_pair(parent2, parent1);

            /* The range that will go from parent1 -> child1. (Not using the constructor is intentional.) */
            unordered_set<size_t> range1;
            for (auto gene = parent1.chromosome.begin() + idx1; gene != parent1.chromosome.begin() + idx2 + 1; gene++) range1.insert(*gene);
            /* The range that will go from parent2 -> child2. */
            unordered_set<size_t> range2;
            for (auto gene = parent2.chromosome.begin() + idx1; gene != parent2.chromosome.begin() + idx2 + 1; gene++) range2.insert(*gene);

            /* Gather genes not in the range from the other parent. */
            vector<size_t> seg1;	/* Segment gathered from parent2 -> child1. */
            seg1.reserve(parent2.chromosome.size() - idx2 + idx1 - 1);
            vector<size_t> seg2;	/* Segment gathered from parent1 -> child2. */
            seg2.reserve(parent1.chromosome.size() - idx2 + idx1 - 1);
            for (size_t i = 0; i < parent1.chromosome.size(); i++)
            {
                /* If this gene of parent2 is not in the range from parent1, add it to seg1. */
                if (!range1.contains(parent2.chromosome[i])) seg1.push_back(parent2.chromosome[i]);
                /* If this gene of parent1 is not in the range from parent2, add it to seg2. */
                if (!range2.contains(parent1.chromosome[i])) seg2.push_back(parent1.chromosome[i]);
            }

            /* Construct the children. child1 = seg1 + p1_range, child2 = seg2 + p2_range */
            child1.chromosome.insert(child1.chromosome.end(), seg1.begin(), seg1.begin() + idx1);
            child1.chromosome.insert(child1.chromosome.end(), parent1.chromosome.begin() + idx1, parent1.chromosome.begin() + idx2 + 1);
            child1.chromosome.insert(child1.chromosome.end(), seg1.begin() + idx1, seg1.end());

            child2.chromosome.insert(child2.chromosome.end(), seg2.begin(), seg2.begin() + idx1);
            child2.chromosome.insert(child2.chromosome.end(), parent2.chromosome.begin() + idx1, parent2.chromosome.begin() + idx2 + 1);
            child2.chromosome.insert(child2.chromosome.end(), seg2.begin() + idx1, seg2.end());
        }
        else  /* No crossover. */
        {
            child1 = parent1;
            child2 = parent2;
        }

        return make_pair(child1, child2);
    }

    inline PermutationGA::CandidatePair PermutationGA::cycleCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
    {
        using namespace std;
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(0.0 <= pc && pc <= 1.0);

        Candidate child1(parent1), child2(parent2);

        /* Crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            /* Identify all cycles. */
            vector<vector<size_t>> cycles;
            /* Copies of parent chromosomes so they can be changed without issues. */
            vector<size_t> chrom1 = parent1.chromosome;
            vector<size_t> chrom2 = parent2.chromosome;
            while (!chrom1.empty())
            {
                /* Identify 1 cycle. */
                vector<size_t> cycle;
                /* Always start the cycle at chrom1[0]. */
                size_t pos = 0;
                size_t top_val = chrom1[pos];
                cycle.push_back(top_val);

                while (chrom2[pos] != chrom1[0])
                {
                    /* Look for the bottom value at this pos (chrom2[pos]) in chrom1. This is the new pos and top value. */
                    pos = static_cast<size_t>(find(chrom1.begin(), chrom1.end(), chrom2[pos]) - chrom1.begin());
                    top_val = chrom1[pos];
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
                /* Odd cycles were already handled when initializing the children. */
            }
            child1.is_evaluated = false;
            child2.is_evaluated = false;
        }

        return make_pair(child1, child2);
    }

    inline PermutationGA::CandidatePair PermutationGA::edgeCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
    {
        using namespace std;
        using NList = vector<unordered_set<size_t>>;

        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(0.0 <= pc && pc <= 1.0);

        Candidate child1, child2;
        child1.chromosome.reserve(parent1.chromosome.size());
        child2.chromosome.reserve(parent2.chromosome.size());

        /* Crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
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
            NList nl2 = nl1;	/* Copy for child2. */

            /* Generate child1. */
            size_t X = parent1.chromosome[0];
            vector<size_t> not_in_child = parent1.chromosome;
            while (child1.chromosome.size() != parent1.chromosome.size())
            {
                /* Append X to the child, and remove X from all neighbour lists. */
                child1.chromosome.push_back(X);
                not_in_child.erase(remove(not_in_child.begin(), not_in_child.end(), X), not_in_child.end());
                for (auto& elem : nl1)
                {
                    elem.erase(X);
                }

                /* Determine next X that will be added to the child. */
                if (child1.chromosome.size() != parent1.chromosome.size())
                {
                    /* If X's neighbour list is empty, X = random node not already in child. */
                    if (nl1[X].empty())
                    {
                        X = not_in_child[rng::randomIdx(not_in_child.size())];
                    }
                    else /* X's neighbour list is not empty, X = neighbour of X with fewest neighbours (random if tie). */
                    {
                        /* Find X's neighbour with fewest neighbours. */
                        size_t nb = *min_element(nl1[X].begin(), nl1[X].end(),
                        [&nl1](const size_t& lhs, const size_t& rhs)
                        {
                            return (nl1[lhs].size() < nl1[rhs].size());
                        });
                        size_t min_neighbour_count = nl1[nb].size();

                        /* Determine possible nodes (neighbours of X with min_neighbour_count num of neighbours). */
                        vector<size_t> possible_nodes;
                        for (const auto& elem : nl1[X])
                        {
                            if (nl1[elem].size() == min_neighbour_count)
                            {
                                possible_nodes.push_back(elem);
                            }
                        }

                        X = possible_nodes[rng::randomIdx(possible_nodes.size())];
                    }
                }
            }

            /* Do same to get child2. */
            X = parent2.chromosome[0];
            not_in_child = parent2.chromosome;
            while (child2.chromosome.size() != parent2.chromosome.size())
            {
                child2.chromosome.push_back(X);
                not_in_child.erase(remove(not_in_child.begin(), not_in_child.end(), X), not_in_child.end());
                for (auto& elem : nl2)
                {
                    elem.erase(X);
                }

                if (child2.chromosome.size() != parent2.chromosome.size())
                {
                    if (nl2[X].empty())
                    {
                        X = not_in_child[rng::randomIdx(not_in_child.size())];
                    }
                    else
                    {
                        size_t nb = *min_element(nl2[X].begin(), nl2[X].end(),
                        [&nl2](const size_t& lhs, const size_t& rhs)
                        {
                            return (nl2[lhs].size() < nl2[rhs].size());
                        });
                        size_t min_neighbour_count = nl2[nb].size();

                        vector<size_t> possible_nodes;
                        for (const auto& elem : nl2[X])
                        {
                            if (nl2[elem].size() == min_neighbour_count)
                            {
                                possible_nodes.push_back(elem);
                            }
                        }
                        X = possible_nodes[rng::randomIdx(possible_nodes.size())];
                    }
                }
            }
        }
        else /* No crossover. */
        {
            child1 = parent1;
            child2 = parent2;
        }

        return make_pair(child1, child2);
    }

    inline PermutationGA::CandidatePair PermutationGA::pmxCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
    {
        using namespace std;
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(0.0 <= pc && pc <= 1.0);

        Candidate child1(parent2), child2(parent1);	/* Init so the last step of the crossover can be skipped. */

        /* Crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            /* Pick a random range of genes. The bounds of the range may be the same, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(parent1.chromosome.size());
            size_t r2 = rng::randomIdx(parent1.chromosome.size());
            const auto [idx1, idx2] = minmax(r1, r2);

            /* Edge case. The entire chromosomes are copied directly. */
            if (idx1 == 0 && idx2 == parent1.chromosome.size() - 1)
            {
                return make_pair(parent1, parent2);
            }

            /* Copy values in the range from the corresponding parent. */
            for (size_t i = idx1; i <= idx2; i++)
            {
                child1.chromosome[i] = parent1.chromosome[i];
                child2.chromosome[i] = parent2.chromosome[i];
            }
            /* Ranges that were copied from parents for fast checking if they contain an element. (Not using the constructor is intentional.) */
            unordered_set<size_t> p1_range;
            for (auto gene = parent1.chromosome.begin() + idx1; gene != parent1.chromosome.begin() + idx2 + 1; gene++)
            {
                p1_range.insert(*gene);
            }
            unordered_set<size_t> p2_range;
            for (auto gene = parent2.chromosome.begin() + idx1; gene != parent2.chromosome.begin() + idx2 + 1; gene++)
            {
                p2_range.insert(*gene);
            }

            /* Get rest of the child genes from the other parents. */
            for (size_t i = idx1; i <= idx2; i++)
            {
                /* Look for genes in parent2 in the same range which haven't been copied to child1 from parent1 (p1_range). */
                if (!p1_range.contains(parent2.chromosome[i]))
                {
                    size_t pos = i;
                    while (idx1 <= pos && pos <= idx2)
                    {
                        /* Look at value in parent1 in this same pos. */
                        size_t val = parent1.chromosome[pos];
                        /* Find this value in parent2. */
                        pos = static_cast<size_t>(find(parent2.chromosome.begin(), parent2.chromosome.end(), val) - parent2.chromosome.begin());
                        /* Keep going until pos is outside the range. */
                    }
                    child1.chromosome[pos] = parent2.chromosome[i];
                }

                /* Same for child2. */
                if (!p2_range.contains(parent1.chromosome[i]))
                {
                    size_t pos = i;
                    while (idx1 <= pos && pos <= idx2)
                    {
                        size_t val = parent2.chromosome[pos];
                        pos = static_cast<size_t>(find(parent1.chromosome.begin(), parent1.chromosome.end(), val) - parent1.chromosome.begin());
                    }
                    child2.chromosome[pos] = parent1.chromosome[i];
                }
            }
            /* Copy any not yet in child positions to the children from the other parents. (Already done at the initialization of the children.) */

            child1.is_evaluated = false;
            child2.is_evaluated = false;
        }

        return make_pair(child1, child2);
    }

    inline void PermutationGA::swapMutate(Candidate& child, double pm)
    {
        assert(0.0 <= pm && pm <= 1.0);

        /* Perform mutation with pm probability. */
        if (rng::randomReal() <= pm)
        {
            /* r1 and r2 might be the same index, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(child.chromosome.size());
            size_t r2 = rng::randomIdx(child.chromosome.size());

            std::swap(child.chromosome[r1], child.chromosome[r2]);

            /* If the indices are different, the child was changed and will need evaluation. */
            if (r1 != r2) child.is_evaluated = false;
        }
    }

    inline void PermutationGA::scrambleMutate(Candidate& child, double pm)
    {
        assert(0.0 <= pm && pm <= 1.0);

        /* Perform mutation with pm probability. */
        if (rng::randomReal() <= pm)
        {
            /* Pick a random range of genes. The bounds may be the same, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(child.chromosome.size());
            size_t r2 = rng::randomIdx(child.chromosome.size());
            auto [idx1, idx2] = std::minmax(r1, r2);

            std::shuffle(child.chromosome.begin() + idx1, child.chromosome.begin() + idx2 + 1, rng::prng);

            /* If the indices are different, the child was very likely changed and will need evaluation. */
            if (r1 != r2) child.is_evaluated = false;
        }
    }

    inline void PermutationGA::inversionMutate(Candidate& child, double pm)
    {
        assert(0.0 <= pm && pm <= 1.0);

        /* Perform mutation with pm probability. */
        if (rng::randomReal() <= pm)
        {
            /* Pick a random range of genes. The bounds of the range may be the same, but its rare for long chromosomes. */
            size_t r1 = rng::randomIdx(child.chromosome.size());
            size_t r2 = rng::randomIdx(child.chromosome.size());
            auto [idx1, idx2] = std::minmax(r1, r2);

            std::reverse(child.chromosome.begin() + idx1, child.chromosome.begin() + idx2 + 1);

            /* If the indices are different, the child was changed and will need evaluation. */
            if (r1 != r2) child.is_evaluated = false;
        }
    }

} // namespace genetic_algorithm


#endif // !GA_PERMUTATION_GA_H