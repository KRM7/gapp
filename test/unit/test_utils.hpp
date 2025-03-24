/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_UTILS_HPP
#define GA_TEST_UTILS_HPP

#include "core/fitness_function.hpp"
#include "core/candidate.hpp"
#include <vector>
#include <cstddef>

using namespace gapp;

template<typename T>
class DummyFitnessFunction final : public FitnessFunctionBase<T>
{
public:
    using Type = FitnessFunctionInfo::Type;

    explicit DummyFitnessFunction(size_t chrom_len, size_t nobj = 1, Type type = Type::Static) :
        FitnessFunctionBase<T>(chrom_len, type), nobj_(nobj)
    {}

    template<size_t N>
    explicit DummyFitnessFunction(std::array<size_t, N> chrom_lens, size_t nobj = 1, Type type = Type::Static) :
        FitnessFunctionBase<T>(chrom_lens, type), nobj_(nobj)
    {}

private:
    size_t nobj_;

    FitnessVector invoke(const Candidate<T>&) const override
    {
        return FitnessVector(nobj_, 0.0); // NOLINT(*return-braced-init-list)
    }
};

#endif // !GA_TEST_UTILS_HPP
