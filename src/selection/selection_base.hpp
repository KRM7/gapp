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

#include <vector>
#include <cstddef>

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

        virtual void prepare(const GA<T>& ga) = 0;
        virtual Candidate<T> select(const GA<T>& ga) = 0;

    private:
        std::vector<double> pdf_;   /* Discrete selection probability distribution function of the population. */
        std::vector<double> cdf_;   /* Discrete cumulative distribution function of the population. */
    };

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_BASE_HPP