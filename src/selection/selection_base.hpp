#ifndef GA_SELECTION_BASE_HPP
#define GA_SELECTION_BASE_HPP

#include "../base_ga.decl.hpp"
#include "../concepts.hpp"

#include <vector>

namespace genetic_algorithm::selection
{
    template<gene T>
    class Selection
    {
    public:
        Selection() = default;
        virtual ~Selection() = default;

        virtual void prepare(const GA<T>& ga) = 0;
        virtual void select(const GA<T>& ga) = 0;

    private:
        std::vector<double> pdf_;   /* Discrete selection probability distribution function of the population. */
        std::vector<double> cdf_;   /* Discrete cumulative distribution function of the population. */
    };

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_BASE_HPP