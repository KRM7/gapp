/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_UTILS_HPP
#define GA_TEST_UTILS_HPP

#include "core/fitness_function.hpp"
#include "population/candidate.hpp"
#include <vector>
#include <cstddef>

using namespace gapp;

template<typename T>
class DummyFitnessFunction final : public FitnessFunctionBase<T>
{
public:
    explicit DummyFitnessFunction(size_t chrom_len, size_t nobj = 1, bool var_len = false, bool dynamic = false) :
        FitnessFunctionBase<T>(chrom_len, var_len, dynamic), nobj_(nobj) {}
private:
    FitnessVector invoke(const Chromosome<T>&) const override { return FitnessVector(nobj_, 0.0); }
    size_t nobj_;
};

#endif // !GA_TEST_UTILS_HPP
