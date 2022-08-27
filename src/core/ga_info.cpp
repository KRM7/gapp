/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "ga_info.hpp"
#include "../utility/utility.hpp"
#include <utility>
#include <stdexcept>

namespace genetic_algorithm
{
    GaInfo::GaInfo(GaInfo&& other) noexcept :
        fitness_matrix_(std::move(other.fitness_matrix_)),
        algorithm_(std::move(other.algorithm_)),
        stop_condition_(std::move(other.stop_condition_)),
        num_fitness_evals_(other.num_fitness_evals_.load()),
        generation_cntr_(other.generation_cntr_),
        num_objectives_(other.num_objectives_),
        chrom_len_(other.chrom_len_),
        population_size_(other.population_size_),
        max_gen_(other.max_gen_),
        dynamic_fitness_(other.dynamic_fitness_),
        variable_chrom_len_(other.variable_chrom_len_),
        keep_all_optimal_sols_(other.keep_all_optimal_sols_),
        can_continue_(other.can_continue_)
    {}

    GaInfo& GaInfo::operator=(GaInfo&& other) noexcept
    {
        if (&other != this)
        {
            fitness_matrix_ = std::move(other.fitness_matrix_);

            algorithm_ = std::move(other.algorithm_);
            stop_condition_ = std::move(other.stop_condition_);

            num_fitness_evals_ = other.num_fitness_evals_.load();
            generation_cntr_ = other.generation_cntr_;
            num_objectives_ = other.num_objectives_;

            chrom_len_ = other.chrom_len_;
            population_size_ = other.population_size_;
            max_gen_ = other.max_gen_;

            dynamic_fitness_ = other.dynamic_fitness_;
            variable_chrom_len_ = other.variable_chrom_len_;
            keep_all_optimal_sols_ = other.keep_all_optimal_sols_;
            can_continue_ = other.can_continue_;
        }

        return *this;
    }

    GaInfo::~GaInfo() = default;

    GaInfo::GaInfo(size_t chrom_len)
        : GaInfo(DEFAULT_POPSIZE, chrom_len)
    {}

    GaInfo::GaInfo(size_t population_size, size_t chrom_len)
    {
        this->population_size(population_size);
        this->chrom_len(chrom_len);
    }

    void GaInfo::chrom_len(size_t len)
    {
        if (len == 0) GA_THROW(std::invalid_argument, "The chromosome length must be at least 1.");

        chrom_len_ = len;
    }

    void GaInfo::population_size(size_t size)
    {
        if (size == 0) GA_THROW(std::invalid_argument, "The population size must be at least 1.");

        population_size_ = size;
    }

    void GaInfo::max_gen(size_t max_gen)
    {
        if (max_gen == 0) GA_THROW(std::invalid_argument, "The number of generations must be at least 1.");

        max_gen_ = max_gen;
    }

    void GaInfo::num_objectives(size_t n)
    {
        if (n == 0) GA_THROW(std::invalid_argument, "There must be at least 1 objective function.");

        num_objectives_ = n;
    }

    void GaInfo::stop_condition(StopConditionFunction f)
    {
        if (!f) GA_THROW(std::invalid_argument, "The stop condition function can't be a nullptr.");

        stop_condition_ = std::make_unique<stopping::dtl::Lambda>(std::move(f));
    }

} // namespace genetic_algorithm