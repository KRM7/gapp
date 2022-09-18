/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_CONE_TREE_HPP
#define GA_UTILITY_CONE_TREE_HPP

#include <vector>
#include <type_traits>
#include <cstddef>
#include "math.hpp"
#include "matrix.hpp"

namespace genetic_algorithm::detail
{
    using math::Point;

    /* 
    * This data structure is used to solve the maximum inner product search problem in the NSGA-III algorithm
    * when looking for the nearest reference point to a candidate solution.
    *  
    * The implementation is based on:
    *  Ram, Parikshit, and Alexander G. Gray. "Maximum inner-product search using cone trees.", 2012.
    */
    class ConeTree
    {
    public:
        using iterator       = Matrix<double>::iterator;
        using const_iterator = Matrix<double>::const_iterator;

        struct FindResult
        {
            const_iterator elem;
            double prod;
        };

        struct Node
        {
            Point center;
            double radius;
            size_t first;
            size_t last;

            size_t left;
            size_t right;
        };

        ConeTree() = default;

        template<std::input_iterator Iter>
        requires std::is_same_v<std::remove_cvref_t<typename std::iterator_traits<Iter>::value_type>, Point>
        ConeTree(Iter first, Iter last);

        size_t size() const noexcept { return points_.size(); }

        Matrix<double>& data() noexcept { return points_; }
        const Matrix<double>& data() const noexcept { return points_; }

        /* Returns the closest point in the tree to the point, and its distance. */
        FindResult findBestMatch(const Point& query_point) const;

    private:
        Matrix<double> points_;
        std::vector<Node> nodes_;

        /* The maximum number of elements allowed in a leaf node. */
        inline static constexpr size_t MAX_LEAF_ELEMENTS = 22;

        /* Expand the tree from the root. */
        void buildTree();
    };

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

namespace genetic_algorithm::detail
{
    template<std::input_iterator Iter>
    requires std::is_same_v<std::remove_cvref_t<typename std::iterator_traits<Iter>::value_type>, Point>
    ConeTree::ConeTree(Iter first, Iter last) :
        points_(0, (*first).size())
    {
        points_.reserve(last - first, points_.ncols());
        while (first != last) points_.push_back(*first++);

        nodes_.reserve(4 * points_.size() / MAX_LEAF_ELEMENTS);
        Node root{ .first = 0, .last = points_.size() };
        nodes_.push_back(root);

        buildTree();
    }
}

#endif // !GA_UTILITY_CONE_TREE_HPP