/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "ga_info.hpp"
#include "../algorithm/algorithm_base.hpp"
#include "../stop_condition/stop_condition.hpp"
#include "../utility/utility.hpp"
#include <atomic>
#include <memory>
#include <utility>

namespace gapp
{
    GaInfo::GaInfo(GaInfo&&) noexcept            = default;
    GaInfo& GaInfo::operator=(GaInfo&&) noexcept = default;

    GaInfo::~GaInfo()                            = default;


    GaInfo::GaInfo(Positive<size_t> population_size, std::unique_ptr<algorithm::Algorithm> algorithm, std::unique_ptr<stopping::StopCondition> stop_condition) noexcept :
        algorithm_(std::move(algorithm)), stop_condition_(std::move(stop_condition)), population_size_(population_size)
    {
        GAPP_ASSERT(stop_condition_, "The stop condition can't be a nullptr.");
    }

    size_t GaInfo::num_fitness_evals() const noexcept
    {
        std::atomic_ref num_fitness_evals{ num_fitness_evals_ };
        return num_fitness_evals.load(std::memory_order_acquire);
    }

    void GaInfo::algorithm(std::unique_ptr<algorithm::Algorithm> f)
    {
        GAPP_ASSERT(f, "The algorithm can't be a nullptr.");

        algorithm_ = std::move(f);
        use_default_algorithm_ = false;
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