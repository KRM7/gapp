/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ALGORITHM_BASE_FWD_HPP
#define GA_ALGORITHM_ALGORITHM_BASE_FWD_HPP

#include <concepts>

namespace genetic_algorithm::algorithm
{
    class Algorithm;

    /** Algorithm types. */
    template<typename T>
    concept AlgorithmType = requires
    {
        requires std::derived_from<T, Algorithm>;
    };

} // namespace genetic_algorithm::algorithm

#endif //!GA_ALGORITHM_ALGORITHM_BASE_FWD_HPP