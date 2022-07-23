/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_SELECTION_BASE_FWD_HPP
#define GA_SELECTION_BASE_FWD_HPP

#include <concepts>

namespace genetic_algorithm::selection
{
    class Selection;

    template<typename T>
    concept SelectionMethod = requires
    {
        requires std::derived_from<T, Selection>;
    };

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_BASE_FWD_HPP