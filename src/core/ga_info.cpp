/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "ga_info.hpp"
#include "../algorithm/single_objective.hpp"
#include "../stop_condition/stop_condition.hpp"
#include "../utility/utility.hpp"
#include <atomic>
#include <memory>
#include <utility>

namespace gapp
{
    const small_vector<size_t> empty;

    GaInfo::GaInfo(GaInfo&&) noexcept            = default;
    GaInfo& GaInfo::operator=(GaInfo&&) noexcept = default;

    GaInfo::~GaInfo() noexcept                   = default;


    GaInfo::GaInfo(Positive<size_t> population_size, std::unique_ptr<algorithm::Algorithm> algorithm, std::unique_ptr<stopping::StopCondition> stop_condition) :
        algorithm_(std::move(algorithm)), stop_condition_(std::move(stop_condition)), population_size_(population_size), use_default_algorithm_(!algorithm_)
    {
        if (!algorithm_) algorithm_ = std::make_unique<algorithm::SingleObjective>();
        if (!stop_condition_) stop_condition_ = std::make_unique<stopping::NoEarlyStop>();
    }

    const FitnessFunctionInfo* GaInfo::fitness_function() const& noexcept
    {
        return fitness_function_.get();
    }

    const small_vector<size_t>& GaInfo::chrom_lens() const noexcept
    {
        return fitness_function_ ? fitness_function_->chrom_lens() : empty;
    }

    size_t GaInfo::num_fitness_evals() const noexcept
    {
        return std::atomic_ref{ num_fitness_evals_ }.load(std::memory_order_acquire);
    }

    void GaInfo::algorithm(std::unique_ptr<algorithm::Algorithm> f)
    {
        use_default_algorithm_ = !f;
        algorithm_ = f ? std::move(f) : std::make_unique<algorithm::SingleObjective>();
    }

    void GaInfo::algorithm(std::nullptr_t)
    {
        use_default_algorithm_ = true;
        algorithm_ = std::make_unique<algorithm::SingleObjective>();
    }

    void GaInfo::stop_condition(std::unique_ptr<stopping::StopCondition> f)
    {
        stop_condition_ = f ? std::move(f) : std::make_unique<stopping::NoEarlyStop>();
    }

    void GaInfo::stop_condition(std::nullptr_t)
    {
        stop_condition_ = std::make_unique<stopping::NoEarlyStop>();
    }

    void GaInfo::stop_condition(StopConditionCallable f)
    {
        stop_condition_ = std::make_unique<stopping::Lambda>(std::move(f));
    }

} // namespace gapp
