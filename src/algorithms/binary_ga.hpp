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
    class BinaryGA : public GA<char>
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
        Candidate generateCandidate() const override;
    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include "../rng.hpp"
#include <cassert>
#include <cstdlib>

namespace genetic_algorithm
{
    inline BinaryGA::BinaryGA(size_t chrom_len, FitnessFunction fitness_function)
        : GA(chrom_len, fitness_function) 
    {
    }

    inline BinaryGA::Candidate BinaryGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(char(rng::randomBool()));
        }

        return sol;
    }

} // namespace genetic_algorithm

#endif // !GA_BINARY_GA_HPP