/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_GENETIC_ALGORITHM_HPP
#define GA_GENETIC_ALGORITHM_HPP

#include "utility/rng.hpp"
#include "algorithm/algorithm.hpp"
#include "crossover/crossover.hpp"
#include "mutation/mutation.hpp"
#include "stop_condition/stop_condition.hpp"
#include "stop_condition/composite.hpp"
#include "core/ga_base.hpp"
#include "encoding/binary.hpp"
#include "encoding/real.hpp"
#include "encoding/permutation.hpp"
#include "encoding/integer.hpp"

namespace ga = genetic_algorithm;

#endif // !GA_GENETIC_ALGORITHM_HPP