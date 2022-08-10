/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_POP_UPDATE_HPP
#define GA_ALGORITHM_POP_UPDATE_HPP

#include "../population/population.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;
}

/** Methods used to update the populations in each generation of the algorithms. */
namespace genetic_algorithm::update
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /**
    * A population update method that selects only the child Candidates from the
    * combined parent and child populations, and uses these as the population of
    * the next generation of the algorithm. \n
    * If the number of children is greater than the population size used in the algorithm,
    * only the first pop_size children will be selected.
    */
    class KeepChildren final
    {
    public:
        /**
        * Take a FitnessMatrix of the combined parent and child solutions specified by
        * the range [first, last), and return the indices of the fitness vectors that
        * belong to the children. \n
        * 
        * In the range [first, last), the [first, children_first) range should be the fitness matrix
        * of the parent population, while the range [children_first, last) should be the fitness matrix
        * of the child population. \n
        * 
        * For the returned indices, it is assumed that the index of @p first is 0.
        * 
        * @param ga The genetic algorithm this operator is used in.
        * @param first Iterator pointing to the fitness vector of the first parent.
        * @param children_first Iterator pointing to the fitness vector of the first child.
        * @param last Iterator pointing to the last element of the fitness matrix.
        * @returns The indices of the Candidates that will make up the next generation of the algorithm.
        */
        std::vector<size_t> operator()(const GaInfo& ga,
                                       FitnessMatrix::const_iterator first,
                                       FitnessMatrix::const_iterator children_first,
                                       FitnessMatrix::const_iterator last);
    };

    /**
    * A population update method that selects the candidates of the next generation
    * using elitism. \n
    * Of the combined parent and child populations, the N best candidates of the parent
    * population are carried over to the next population, while the remaining
    * (pop_size - N) slots are filled by the child solutions. \n
    * If N = 0, this is equivalent to only keeping the children for the next
    * generation (KeepChildren).
    */
    class Elitism final
    {
    public:
        /**
        * Create an elitist population update operator.
        * 
        * @param n The number of the best parents that will be carried over to the next
        *          generation of the algorithm.
        */
        Elitism(size_t n = 1);

        /**
        * Set the number of elites used to @p n.
        * 
        * @param n The number of the best parents that will be carried over to the next
        *          generation of the algorithm.
        */
        void elite_num(size_t n) noexcept;

        /** @returns The number of elites used. */
        [[nodiscard]]
        size_t elite_num() noexcept { return n_; }

        /**
        * Take a FitnessMatrix of the combined parent and child solutions specified by
        * the range [first, last), and return the indices of the fitness vectors that
        * are going to make up the population of the next generation of the algorithm. \n
        *
        * In the range [first, last), the [first, children_first) range should be the fitness matrix
        * of the parent population, while the range [children_first, last) should be the fitness matrix
        * of the child population. \n
        *
        * For the returned indices, it is assumed that the index of @p first is 0.
        *
        * @param ga The genetic algorithm this operator is used in.
        * @param first Iterator pointing to the fitness vector of the first parent.
        * @param children_first Iterator pointing to the fitness vector of the first child.
        * @param last Iterator pointing to the last element of the fitness matrix.
        * @returns The indices of the Candidates that will make up the next generation of the algorithm.
        */
        std::vector<size_t> operator()(const GaInfo& ga,
                                       FitnessMatrix::const_iterator parents_first,
                                       FitnessMatrix::const_iterator children_first,
                                       FitnessMatrix::const_iterator children_last);

    private:
        size_t n_;
    };

    /**
    * A population update method that selects the best (pop_size) candidates of the combined
    * parent and child populations (assuming fitness maximization), and uses these as the
    * candidates of the next generation of the algorithm.
    */
    class KeepBest final
    {
    public:
        /**
        * Take a FitnessMatrix of the combined parent and child solutions specified by
        * the range [first, last), and return the indices of the fitness vectors that
        * are going to make up the population of the next generation of the algorithm. \n
        *
        * In the range [first, last), the [first, children_first) range should be the fitness matrix
        * of the parent population, while the range [children_first, last) should be the fitness matrix
        * of the child population. \n
        *
        * For the returned indices, it is assumed that the index of @p first is 0.
        *
        * @param ga The genetic algorithm this operator is used in.
        * @param first Iterator pointing to the fitness vector of the first parent.
        * @param children_first Iterator pointing to the fitness vector of the first child.
        * @param last Iterator pointing to the last element of the fitness matrix.
        * @returns The indices of the Candidates that will make up the next generation of the algorithm.
        */
        std::vector<size_t> operator()(const GaInfo& ga,
                                       FitnessMatrix::const_iterator parents_first,
                                       FitnessMatrix::const_iterator children_first,
                                       FitnessMatrix::const_iterator children_last);
    };

} // namespace genetic_algorithm::update


#endif // !GA_ALGORITHM_POP_UPDATE_HPP