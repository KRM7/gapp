/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_NEIGHBOUR_LIST_HPP
#define GAPP_CROSSOVER_NEIGHBOUR_LIST_HPP

#include "../core/candidate.hpp"
#include "../utility/iterators.hpp"
#include "../utility/utility.hpp"
#include <array>
#include <cstddef>

namespace gapp::crossover::impl
{
    class NeighbourList : detail::iterator_interface<NeighbourList>
    {
    public:
        constexpr auto begin() noexcept { return neighbours_.begin(); }
        constexpr auto end() noexcept { return neighbours_.end(); }
        constexpr auto begin() const noexcept { return neighbours_.begin(); }
        constexpr auto end() const noexcept { return neighbours_.end(); }

        constexpr size_t size() const noexcept { return size_; }
        constexpr bool empty() const noexcept { return size() == 0; }

        void insert(PermutationGene value)
        {
            GAPP_ASSERT(value != EMPTY);

            /* Assume that EMPTY values are at the back. */
            for (PermutationGene& neighbour : neighbours_)
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

        void remove(PermutationGene value)
        {
            GAPP_ASSERT(value != EMPTY);

            for (PermutationGene& neighbour : neighbours_)
            {
                if (neighbour == value)
                {
                    neighbour = EMPTY;
                    size_--;
                }
            }
        }

        static constexpr PermutationGene EMPTY = PermutationGene(-1);

    private:
        std::array<PermutationGene, 4> neighbours_{ EMPTY, EMPTY, EMPTY, EMPTY };
        size_t size_ = 0;
    };

    using NeighbourLists = std::vector<NeighbourList>;

    inline NeighbourLists makeNeighbourLists(const Chromosome<PermutationGene>& chrom1, const Chromosome<PermutationGene>& chrom2)
    {
        GAPP_ASSERT(chrom1.size() == chrom2.size());

        NeighbourLists nb_lists(chrom1.size());

        nb_lists[chrom1.front()].insert(chrom1[1]);
        nb_lists[chrom2.front()].insert(chrom2[1]);

        for (size_t i = 1; i < chrom1.size() - 1; i++)
        {
            nb_lists[chrom1[i]].insert(chrom1[i - 1]);
            nb_lists[chrom1[i]].insert(chrom1[i + 1]);

            nb_lists[chrom2[i]].insert(chrom2[i - 1]);
            nb_lists[chrom2[i]].insert(chrom2[i + 1]);
        }

        nb_lists[chrom1.back()].insert(*(chrom1.end() - 2));
        nb_lists[chrom2.back()].insert(*(chrom2.end() - 2));

        return nb_lists;
    }

} // namespace gapp::crossover::impl

#endif // !GAPP_CROSSOVER_NEIGHBOUR_LIST_HPP
