/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_UTILS_HPP
#define GA_TEST_UTILS_HPP

#include "core/fitness_function.hpp"
#include "population/candidate.hpp"
#include <vector>
#include <cstddef>

using namespace genetic_algorithm;

template<typename T>
class DummyFitnessFunction final : public FitnessFunction<T>
{
public:
    explicit DummyFitnessFunction(size_t chrom_len, size_t nobj = 1, bool var_len = false, bool dynamic = false) :
        FitnessFunction<T>(chrom_len, nobj, var_len, dynamic) {}
private:
    detail::FitnessVector invoke(const Chromosome<T>&) const override { return std::vector(this->num_objectives(), 0.0); }
};

#endif // !GA_TEST_UTILS_HPP
