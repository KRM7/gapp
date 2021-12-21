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

#ifndef GA_MUTATION_INTEGER_HPP
#define GA_MUTATION_INTEGER_HPP

#include "mutation_base.hpp"

#include <cstddef>

namespace genetic_algorithm::mutation::integer
{
    class Uniform : public Mutation<size_t>
    {
    public:
        Uniform(size_t base, double pc = 0.01);

        void base(size_t base);

        [[nodiscard]]
        size_t base() const noexcept { return base_; }
    private:
        size_t base_;

        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::integer

#endif // !GA_MUTATION_INTEGER_HPP