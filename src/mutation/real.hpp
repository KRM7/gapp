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

#ifndef GA_MUTATION_REAL_HPP
#define GA_MUTATION_REAL_HPP

#include "mutation_base.hpp"

#include <vector>
#include <utility>

namespace genetic_algorithm::mutation::real
{
    class Uniform final : public BoundedMutation<double>
    {
    public:
        using BoundedMutation::BoundedMutation;
    private:
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;
    };

    class NonUniform final : public BoundedMutation<double>
    {
    public:
        NonUniform(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double beta = 2.0);
        void beta(double b);
        [[nodiscard]] double beta() const noexcept { return beta_; }
    private:
        double beta_ = 2.0;
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;
    };

    class Gauss final : public BoundedMutation<double>
    {
    public:
        Gauss(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double sigma = 6.0);
        void sigma(double sigma);
        [[nodiscard]] double sigma() const noexcept { return sigma_; }
    private:
        double sigma_ = 6.0;
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;
    };

    class Polynomial final : public BoundedMutation<double>
    {
    public:
        Polynomial(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double param = 40.0);
        void param(double param);
        [[nodiscard]] double param() const noexcept { return param_; }
    private:
        double param_ = 40.0;
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;
    };

    class Boundary final : public BoundedMutation<double>
    {
    public:
        using BoundedMutation::BoundedMutation;
    private:
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::real

#endif // !GA_MUTATION_REAL_HPP