/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

#ifndef GA_SELECTION_BASE_HPP
#define GA_SELECTION_BASE_HPP

#include "../concepts.hpp"
#include "../population.hpp"

namespace genetic_algorithm
{
    template<typename T>
    class GA;

} // namespace genetic_algorithm

namespace genetic_algorithm::selection
{
    template<gene T>
    class Selection
    {
    public:
        Selection() = default;
        Selection(const GA<T>& ga) {};
        virtual ~Selection() = default;

        virtual void prepare(const GA<T>& ga, const Population<T>& pop) = 0;
        virtual Candidate<T> select(const GA<T>& ga, const Population<T>& pop) = 0;
        virtual Population<T> nextPopulation(const GA<T>& ga, Population<T>& old_pop, CandidateVec<T>& children) const;
    };

} // namespace genetic_algorithm::selection


/* IMPLEMENTATION */

#include <algorithm>
#include <iterator>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::selection
{
    template<gene T>
    inline Population<T> Selection<T>::nextPopulation(const GA<T>&, Population<T>& old_pop, CandidateVec<T>& children) const
    {
        assert(!children.empty());

        size_t pop_size = old_pop.size();
        old_pop.insert(old_pop.end(), std::make_move_iterator(children.begin()),
                                      std::make_move_iterator(children.end()));

        assert(std::all_of(old_pop.begin(), old_pop.end(), [](const Candidate<T>& sol) { return sol.is_evaluated && sol.fitness.size() > 0; }));

        std::partial_sort(old_pop.begin(), old_pop.begin() + pop_size, old_pop.end(),
        [](const Candidate<T>& lhs, const Candidate<T>& rhs)
        {
            return detail::paretoCompareLess(rhs.fitness, lhs.fitness);
        });
        old_pop.resize(pop_size);

        return old_pop;
    }

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_BASE_HPP