/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_NEIGHBOUR_LIST_HPP
#define GAPP_CROSSOVER_NEIGHBOUR_LIST_HPP

#include "../core/candidate.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/iterators.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/utility.hpp"
#include <array>
#include <vector>
#include <unordered_map>
#include <type_traits>
#include <cstddef>

namespace gapp::crossover::dtl
{
    template<typename T>
    class NeighbourList : public detail::container_interface<NeighbourList<T>>
    {
    public:
        using value_type             = T;
        using reference              = T&;
        using const_reference        = const T&;

        using size_type              = std::size_t;
        using difference_type        = std::ptrdiff_t;

        using iterator               = typename small_vector<value_type>::iterator;
        using const_iterator         = typename small_vector<value_type>::const_iterator;
        using reverse_iterator       = typename small_vector<value_type>::reverse_iterator;
        using const_reverse_iterator = typename small_vector<value_type>::const_reverse_iterator;

        constexpr auto begin() noexcept { return neighbours_.begin(); }
        constexpr auto end() noexcept { return neighbours_.end(); }
        constexpr auto begin() const noexcept { return neighbours_.begin(); }
        constexpr auto end() const noexcept { return neighbours_.end(); }

        constexpr void add(const T& value)
        {
            if (detail::contains(neighbours_.begin(), neighbours_.end(), value)) return;
            neighbours_.push_back(value);
        }

        constexpr void remove(const T& value)
        {
            detail::erase_first_stable(neighbours_, value);
        }

    private:
        small_vector<T, 4> neighbours_;
    };

    template<std::unsigned_integral T>
    class NeighbourList<T> : public detail::iterator_interface<NeighbourList<T>>
    {
    public:
        constexpr auto begin() noexcept { return neighbours_.begin(); }
        constexpr auto end() noexcept { return neighbours_.end(); }
        constexpr auto begin() const noexcept { return neighbours_.begin(); }
        constexpr auto end() const noexcept { return neighbours_.end(); }

        constexpr size_t size() const noexcept { return size_; }
        constexpr bool empty() const noexcept { return size() == 0; }

        void add(T value)
        {
            GAPP_ASSERT(value != EMPTY);

            /* Assume that EMPTY values are at the back. */
            for (T& neighbour : neighbours_)
            {
                if (neighbour == value) return;
                if (neighbour == EMPTY)
                {
                    neighbour = value;
                    size_++;
                    return;
                }
            }
            GAPP_UNREACHABLE();
        }

        void remove(T value)
        {
            GAPP_ASSERT(value != EMPTY);

            for (T& neighbour : neighbours_)
            {
                if (neighbour == value)
                {
                    neighbour = EMPTY;
                    size_--;
                }
            }
        }

        static constexpr T EMPTY = T(-1);

    private:
        std::array<T, 4> neighbours_{ EMPTY, EMPTY, EMPTY, EMPTY };
        size_t size_ = 0;
    };


    template<typename T>
    using NeighbourLists =
        std::conditional_t<std::is_unsigned_v<T>, std::vector<NeighbourList<T>>, std::unordered_map<T, NeighbourList<T>>>;


    template<typename T>
    NeighbourLists<T> makeNeighbourLists(const Chromosome<T>& chrom1, const Chromosome<T>& chrom2)
    {
        GAPP_ASSERT(chrom1.size() == chrom2.size());

        NeighbourLists<T> nb_lists(chrom1.size());

        nb_lists[chrom1.front()].add(chrom1[1]);
        nb_lists[chrom2.front()].add(chrom2[1]);

        for (size_t i = 1; i < chrom1.size() - 1; i++)
        {
            nb_lists[chrom1[i]].add(chrom1[i - 1]);
            nb_lists[chrom1[i]].add(chrom1[i + 1]);

            nb_lists[chrom2[i]].add(chrom2[i - 1]);
            nb_lists[chrom2[i]].add(chrom2[i + 1]);
        }

        nb_lists[chrom1.back()].add(*(chrom1.end() - 2));
        nb_lists[chrom2.back()].add(*(chrom2.end() - 2));

        return nb_lists;
    }

} // namespace gapp::crossover::dtl

#endif // !GAPP_CROSSOVER_NEIGHBOUR_LIST_HPP
