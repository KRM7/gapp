/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BASE_FWD_HPP
#define GA_CROSSOVER_BASE_FWD_HPP

#include "../population/candidate.hpp"
#include <concepts>

namespace genetic_algorithm::crossover
{
    template<Gene T>
    class Crossover;

    /** Crossover method types. */
    template<typename T, typename G>
    concept CrossoverType = requires
    {
        requires Gene<G>;
        requires std::derived_from<T, Crossover<G>>;
    };

} // namespace genetic_algorithm::crossover

#endif // !GA_CROSSOVER_BASE_FWD_HPP