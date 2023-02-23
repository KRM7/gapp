/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_FWD_HPP
#define GA_MUTATION_BASE_FWD_HPP

#include "../population/candidate.hpp"
#include <concepts>

namespace genetic_algorithm::mutation
{
    template<Gene T>
    class Mutation;

    /** Mutation method types. */
    template<typename T, typename G>
    concept MutationType = requires
    {
        requires Gene<G>;
        requires std::derived_from<T, Mutation<G>>;
    };

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_FWD_HPP