/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ND_SORT_HPP
#define GA_ALGORITHM_ND_SORT_HPP

#include "../core/population.hpp"
#include "../utility/iterators.hpp"
#include <vector>
#include <utility>
#include <cstddef>

namespace gapp::algorithm::dtl
{
    struct FrontElement
    {
        constexpr FrontElement(size_t index, size_t pareto_rank) noexcept :
            idx(index), rank(pareto_rank)
        {}

        friend bool operator==(const FrontElement&, const FrontElement&) = default;

        size_t idx;  // The solution's idx in the fitness matrix
        size_t rank; // The rank of the pareto front the solution belongs to
    };

    class ParetoFrontsRange : public detail::container_interface<ParetoFrontsRange>
    {
    public:
        using value_type      = FrontElement;
        using reference       = FrontElement&;
        using const_reference = const FrontElement&;

        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;

        using iterator               = std::vector<value_type>::iterator;
        using const_iterator         = std::vector<value_type>::const_iterator;
        using reverse_iterator       = std::vector<value_type>::reverse_iterator;
        using const_reverse_iterator = std::vector<value_type>::const_reverse_iterator;

        constexpr ParetoFrontsRange() = default;

        constexpr ParetoFrontsRange(const iterator& first, const iterator& last) noexcept :
            first_(first), last_(last)
        {}

        template<typename Range>
        constexpr explicit ParetoFrontsRange(Range& range) noexcept :
            first_(range.begin()), last_(range.end())
        {}

        constexpr iterator begin() noexcept { return first_; }
        constexpr iterator begin() const noexcept { return first_; }

        constexpr iterator end() noexcept { return last_; }
        constexpr iterator end() const noexcept { return last_; }

    private:
        iterator first_ = {};
        iterator last_  = {};
    };


    class ParetoFronts : public detail::container_interface<ParetoFronts>
    {
    public:
        using value_type      = FrontElement;
        using reference       = FrontElement&;
        using const_reference = const FrontElement&;

        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;

        using iterator               = std::vector<value_type>::iterator;
        using const_iterator         = std::vector<value_type>::const_iterator;
        using reverse_iterator       = std::vector<value_type>::reverse_iterator;
        using const_reverse_iterator = std::vector<value_type>::const_reverse_iterator;

        ParetoFronts(std::vector<FrontElement> fronts);

        iterator begin() noexcept { return elements_.begin(); }
        const_iterator begin() const noexcept { return elements_.begin(); }

        iterator end() noexcept { return elements_.end(); }
        const_iterator end() const noexcept { return elements_.end(); }

        void resize(size_type new_size);

        std::vector<size_t> ranks() const;
        std::vector<ParetoFrontsRange> fronts();
        ParetoFrontsRange partialFront(size_type size);

    private:
        std::vector<FrontElement> elements_;
    };


    std::vector<FrontElement> fastNonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);
    std::vector<FrontElement> dominanceDegreeSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);
    std::vector<FrontElement> efficientNonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    ParetoFronts nonDominatedSort(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);
    ParetoFronts nonDominatedSort(const FitnessMatrix& fmat);

} // namespace gapp::algorithm::dtl

#endif // !GA_ALGORITHM_ND_SORT_HPP