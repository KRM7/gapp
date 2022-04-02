/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "ga_base.decl.hpp"
#include <stdexcept>
#include <cstddef>

namespace genetic_algorithm
{
    GaBase::GaBase(size_t chrom_len)
    {
        if (chrom_len == 0)
        {
            throw std::invalid_argument("The chromosome length must be at least 1.");
        }
        chrom_len_ = chrom_len;
    }

    size_t GaBase::num_fitness_evals() const noexcept
    {
        return num_fitness_evals_.load();
    }

    size_t GaBase::generation_cntr() const noexcept
    {
        return generation_cntr_;
    }

    void GaBase::chrom_len(size_t len)
    {
        if (len == 0)
        {
            throw std::invalid_argument("The chromosome length must be at least 1.");
        }
        chrom_len_ = len;
    }

    size_t GaBase::chrom_len() const noexcept
    {
        return chrom_len_;
    }

    void GaBase::population_size(size_t size)
    {
        if (size == 0)
        {
            throw std::invalid_argument("The population size must be at least 1.");
        }
        population_size_ = size;
    }

    size_t GaBase::population_size() const noexcept
    {
        return population_size_;
    }

    void GaBase::max_gen(size_t max_gen)
    {
        if (max_gen == 0) throw std::invalid_argument("The maximum number of generations must be at least 1.");

        max_gen_ = max_gen;
    }
    
    size_t GaBase::max_gen() const noexcept
    {
        return max_gen_;
    }

    size_t GaBase::num_objectives() const noexcept
    {
        return num_objectives_;
    }

} // namespace genetic_algorithm