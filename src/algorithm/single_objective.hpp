/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_HPP

#include "algorithm_base.hpp"
#include "population_update.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include <cstddef>

namespace genetic_algorithm::algorithm
{
    template<typename Selection, population_update::Updater PopulationUpdate> // concept for selections
    class SingleObjective final : public Algorithm
    {
    public:

        SingleObjective(Selection selection, PopulationUpdate pop_update = population_update::KeepBest{}); // default tourney + keepBest

        void init(const GaInfo& ga) override
        {
            selection_.init(ga); // maybe pass algo as 2nd ctxt, probably unnecessary
        }

        void prepare(const GaInfo& ga, const FitnessMatrix& population_fmat) override
        {
            selection_.prepare(ga, population_fmat); // maybe pass algo as 2nd ctxt
        }

        size_t select(const GaInfo& ga, const FitnessMatrix& population_fmat) override
        {
            return selection_.select(ga, population_fmat); // maybe pass algo as 2nd ctxt
        }

        std::vector<size_t> nextPopulation(const GaInfo& ga, const FitnessMatrix& population_fmat) override // TODO: params should be iterators
        {
            return update_.nextPopulation(ga, population_fmat); // pass comp function (default pareto comp)
        }

        // selection + update getters/setters

    private:
        Selection selection_;
        PopulationUpdate update_;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_SINGLE_OBJECTIVE_HPP