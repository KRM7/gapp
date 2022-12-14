/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_PROBLEMS_TSP_DATA_HPP
#define GA_PROBLEMS_TSP_DATA_HPP

#include <array>

namespace genetic_algorithm::problems
{
    using Coords = std::array<double, 2>;

    constexpr std::array<Coords, 52> tsp52_coords
    {{
        #include "tsp52.txt"
    }};

    constexpr std::array<Coords, 76> tsp76_coords
    {{
        #include "tsp76.txt"
    }};

    constexpr std::array<Coords, 124> tsp124_coords
    {{
        #include "tsp124.txt"
    }};

    constexpr std::array<Coords, 152> tsp152_coords
    {{
        #include "tsp152.txt"
    }};

    constexpr std::array<Coords, 226> tsp226_coords
    {{
        #include "tsp226.txt"
    }};

    constexpr std::array<Coords, 299> tsp299_coords
    {{
        #include "tsp299.txt"
    }};

    constexpr std::array<Coords, 439> tsp439_coords
    {{
        #include "tsp439.txt"
    }};

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_TSP_DATA_HPP