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

#ifndef GA_SELECTION_DECL_HPP
#define GA_SELECTION_DECL_HPP

#include "selection_base.hpp"
#include "../candidate.hpp"
#include "../concepts.hpp"

#include <vector>
#include <utility>
#include <functional>

namespace genetic_algorithm::selection
{
    template<gene T>
    class Roulette final : public Selection<T>
    {
    public:
        using Selection<T>::Selection;

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;
    private:
        std::vector<double> cdf_;
    };

    template<gene T>
    class Tournament final : public Selection<T>
    {
    public:
        using Selection<T>::Selection;

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;
    };

    template<gene T>
    class Rank final : public Selection<T>
    {
    public:
        explicit Rank(double min_weight = 0.1, double max_weight = 1.1);
        explicit Rank(const GA<T>& ga, double min_weight = 0.1, double max_weight = 1.1);

        void min_weight(double min_weight);
        [[nodiscard]]
        double min_weigth() const noexcept { return min_weight_; }

        void max_weight(double max_weight);
        [[nodiscard]]
        double max_weight() const noexcept { return max_weight_; }

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;

    private:
        double min_weight_ = 0.1;
        double max_weight_ = 1.1;

        std::vector<double> cdf_;
    };

    template<gene T>
    class Sigma final : public Selection<T>
    {
    public:
        explicit Sigma(double scale = 3.0);
        explicit Sigma(const GA<T>& ga, double scale = 3.0);

        void scale(double scale);
        [[nodiscard]]
        double scale() const noexcept { return scale_; }

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;

    private:
        double scale_ = 3.0;

        std::vector<double> cdf_;
    };

    template<gene T>
    class Boltzmann final : public Selection<T>
    {
    public:
        explicit Boltzmann(double t_min = 0.25, double t_max = 4.0);
        explicit Boltzmann(const GA<T>& ga, double t_min = 0.25, double t_max = 4.0);

        void t_min(double t_min);
        [[nodiscard]]
        double t_min() const noexcept { return t_min_; }

        void t_max(double t_max);
        [[nodiscard]]
        double t_max() const noexcept { return t_max_; }

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;

    private:
        double t_min_;
        double t_max_;
        //std::function<double(double)> temperature_;

        std::vector<double> cdf_;
    };

    template<gene T>
    class NSGA2 final : public Selection<T>
    {
    public:
        using Selection<T>::Selection;

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;
        Population<T> nextPopulation(const GA<T>& ga, Population<T>& old_pop, CandidateVec<T>& children) const override;

    private:
        std::vector<size_t> rank_;
        std::vector<double> crowd_dist_;
    };

    template<gene T>
    class NSGA3 final : public Selection<T>
    {
    public:
        using Selection<T>::Selection;

        void prepare(const GA<T>& ga, const Population<T>& pop) override;
        Candidate<T> select(const GA<T>& ga, const Population<T>& pop) override;
        Population<T> nextPopulation(const GA<T>& ga, Population<T>& old_pop, CandidateVec<T>& children) const override;

    private:
        std::vector<size_t> rank_;
        std::vector<size_t> ref_idx_;
        std::vector<double> ref_dist_;
        std::vector<size_t> niche_cnt_;
    };

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_DECL_HPP