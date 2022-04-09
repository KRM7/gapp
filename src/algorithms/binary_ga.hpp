/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_BINARY_GA_HPP
#define GA_BINARY_GA_HPP

#include "ga_base.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Standard genetic algorithm with binary encoding. \n
    * (The binary genes are encoded as char types.)
    */
    class BinaryGA final : public GA<char, BinaryGA>
    {
    public:

        /**
        * Basic contructor for the binary GA.
        *
        * @param chrom_len The length of the binary chromosomes.
        * @param fitness_function The fitness function to find the maximum of in the algorithm.
        */
        BinaryGA(size_t chrom_len, FitnessFunction fitness_function);

    private:
        friend class GA<GeneType, BinaryGA>;

        Candidate generateCandidate() const;
    };

} // namespace genetic_algorithm

#endif // !GA_BINARY_GA_HPP