/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_CONE_TREE_HPP
#define GA_UTILITY_CONE_TREE_HPP

#include "math.hpp"
#include "matrix.hpp"
#include <vector>
#include <iterator>
#include <type_traits>
#include <cstddef>

namespace genetic_algorithm::detail
{
    using math::Point;

    /* 
    * This data structure is used for the maximum inner product search in the NSGA-III algorithm,
    * when searching for the nearest reference point to a solution.
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
            Point center  = {};
            double radius = 0.0;
            size_t first  = 0;
            size_t last   = 0;

            size_t left   = 0;
            size_t right  = 0;
        };

        ConeTree() = default;

        template<std::input_iterator Iter>
        requires std::is_same_v<std::remove_cvref_t<typename std::iterator_traits<Iter>::value_type>, Point>
        ConeTree(Iter first, Iter last);


        /* Returns the closest point in the tree to the query point, and its distance. */
        FindResult findBestMatch(const Point& query_point) const;


        const_iterator begin() const noexcept { return points_.begin(); }
        const_iterator end() const noexcept   { return points_.end(); }

        size_t size() const noexcept { return points_.size(); }

        const Matrix<double>& data() const noexcept { return points_; }

    private:
        Matrix<double> points_;
        std::vector<Node> nodes_;

        inline static constexpr size_t MAX_LEAF_ELEMENTS = 22;  /* The maximum number of points in a leaf node. */

        void buildTree();

        iterator node_begin(const Node& node) { return points_.begin() + node.first; }
        iterator node_end(const Node& node) { return points_.begin() + node.last; }

        const_iterator node_begin(const Node& node) const { return points_.begin() + node.first; }
        const_iterator node_end(const Node& node) const { return points_.begin() + node.last; }

        const_iterator node_cbegin(const Node& node) const { return points_.cbegin() + node.first; }
        const_iterator node_cend(const Node& node) const { return points_.cbegin() + node.last; }

        const Node& left_child(const Node& node) const { return nodes_[node.left]; }
        const Node& right_child(const Node& node) const { return nodes_[node.right]; }
    };

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include <cassert>

namespace genetic_algorithm::detail
{
    template<std::input_iterator Iter>
    requires std::is_same_v<std::remove_cvref_t<typename std::iterator_traits<Iter>::value_type>, Point>
    ConeTree::ConeTree(Iter first, Iter last)
    {
        assert(std::distance(first, last) > 0);

        points_ = Matrix(0, (*first).size(), 0.0);
        points_.reserve(last - first, points_.ncols());
        while (first != last) points_.append_row(*first++);

        nodes_.reserve(4 * points_.size() / MAX_LEAF_ELEMENTS);
        Node root{ .first = 0, .last = points_.size() };
        nodes_.push_back(root);

        buildTree();
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_CONE_TREE_HPP