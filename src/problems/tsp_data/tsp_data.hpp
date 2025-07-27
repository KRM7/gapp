/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_PROBLEMS_TSP_DATA_HPP
#define GAPP_PROBLEMS_TSP_DATA_HPP

#include <array>

namespace gapp::problems
{
    using Coords = std::array<double, 2>;

    inline constexpr std::array<Coords, 52> tsp52_coords
    {{
        #include "tsp52.txt"
    }};

    inline  std::array<Coords, 76> tsp76_coords
    {{
        #include "tsp76.txt"
    }};

    inline constexpr std::array<Coords, 124> tsp124_coords
    {{
        #include "tsp124.txt"
    }};

    inline constexpr std::array<Coords, 152> tsp152_coords
    {{
        #include "tsp152.txt"
    }};

    inline constexpr std::array<Coords, 226> tsp226_coords
    {{
        #include "tsp226.txt"
    }};

    inline constexpr std::array<Coords, 299> tsp299_coords
    {{
        #include "tsp299.txt"
    }};

    inline constexpr std::array<Coords, 439> tsp439_coords
    {{
        #include "tsp439.txt"
    }};

} // namespace gapp::problems

#endif // !GAPP_PROBLEMS_TSP_DATA_HPP
