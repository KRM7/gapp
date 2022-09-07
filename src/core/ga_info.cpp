/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "ga_info.hpp"
#include "../algorithm/single_objective.decl.hpp"
#include "../stop_condition/stop_condition_base.hpp"
#include "../algorithm/nsga3.hpp"
#include "../utility/utility.hpp"
#include <utility>
#include <stdexcept>

namespace genetic_algorithm
{
    GaInfo::GaInfo(GaInfo&&) noexcept            = default;
    GaInfo& GaInfo::operator=(GaInfo&&) noexcept = default;
    GaInfo::~GaInfo()                            = default;

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

    void GaInfo::setDefaultAlgorithm()
    {
        num_objectives(findNumObjectives());

        (num_objectives_ == 1) ?
            algorithm(std::make_unique<algorithm::SingleObjective<>>()) :
            algorithm(std::make_unique<algorithm::NSGA3>());
    }

} // namespace genetic_algorithm