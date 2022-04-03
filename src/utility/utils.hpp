/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILS_HPP
#define GA_UTILS_HPP

#include <limits>
#include <execution>
#include <random>

#ifndef GA_EPSILON
#define GA_EPSILON (100 * std::numeric_limits<double>::epsilon())
#endif

#ifndef GA_SEQ_EXECUTION
#define GA_EXECUTION_UNSEQ std::execution::par_unseq
#define GA_EXECUTION_SEQ std::execution::par
#else
#define GA_EXECUTION_UNSEQ std::execution::unseq
#define GA_EXECUTION_SEQ std::execution::seq
#endif

#ifndef GA_PRESET_SEED
#define GA_SEED() (std::random_device{}())
#else
#define GA_SEED() (0x12345678);
#endif

#endif // !GA_UTILS_HPP