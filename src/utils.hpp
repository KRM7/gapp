/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this softwareand associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright noticeand this permission notice shall be included in all
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

/**
* Some utilities.
*
* @file utils.hpp
*/

#ifndef GA_UTILS_HPP
#define GA_UTILS_HPP

#include <cstdlib>
#include <cassert>
#include <limits>
#include <execution>
#include <random>


#define GA_UNREACHABLE (assert(!"Unreachable code."), std::terminate())

#ifndef GA_EPSILON
#define GA_EPSILON (100 * std::numeric_limits<double>::epsilon())
#endif

#ifndef GA_SEQ_EXECUTION
#define GA_EXECUTION_UNSEQ std::execution::par_unseq
#define GA_EXECUTION_SEQ std::execution::par
#else
#define GA_EXECUTION_UNSEQ std::execution::unseq
#define GA_EXECUTION_SEQ std::execution::seq
#endif

#ifndef GA_PRESET_SEED
#define GA_SEED() (std::random_device{}())
#else
#define GA_SEED() (0x12345678);
#endif

#endif // !GA_UTILS_HPP