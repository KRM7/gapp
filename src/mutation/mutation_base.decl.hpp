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

#ifndef GA_MUTATION_BASE_DECL_HPP
#define GA_MUTATION_BASE_DECL_HPP

#include "../candidate.h"

#include <vector>
#include <utility>

namespace genetic_algorithm
{
    template<typename GeneType>
    class GA;

} // namespace genetic_algorithm

namespace genetic_algorithm::mutation
{
    template<regular_hashable GeneType>
    class Mutation
    {
    public:
        explicit Mutation(double pm = 0.01);

        virtual ~Mutation() = default;

        void mutation_rate(double pm);

        [[nodiscard]]
        double mutation_rate() const noexcept { return pm_; }

        void operator()(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const;

    protected:
        virtual void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const = 0;

        double pm_ = 0.01;      /* The crossover rate used with the mutation operator. */
    };

    template<regular_hashable GeneType>
    class BoundedMutation : public Mutation<GeneType>
    {
    public:
        BoundedMutation(const std::vector<std::pair<GeneType, GeneType>>& bounds, double pm = 0.01);

        void bounds(const std::vector<std::pair<GeneType, GeneType>>& bounds);

        [[nodiscard]]
        std::vector<std::pair<GeneType, GeneType>> bounds() const noexcept { return bounds_; }
    protected:
        std::vector<std::pair<GeneType, GeneType>> bounds_;     /* The lower and upper bounds for each gene of the chromosome. */
    };

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_DECL_HPP