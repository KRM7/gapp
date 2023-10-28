/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_CONE_TREE_HPP
#define GA_UTILITY_CONE_TREE_HPP

#include "matrix.hpp"
#include <vector>
#include <span>
#include <type_traits>
#include <cstddef>

namespace gapp::detail
{
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

        using Point = std::vector<double>;

        struct FindResult
        {
            const_iterator elem;
            double prod;
        };

        struct Node
        {
            Point center  = {};
            double radius = 0.0;
            size_t first  = 0;  /* Index of the first point which belongs to the node. */
            size_t last   = 0;  /* Index of the first point which does not belong to the node. */

            size_t left   = 0;  /* Index of the left child node. */
            size_t right  = 0;  /* Index of the right child node. */
        };

        ConeTree() = default;

        explicit ConeTree(std::span<const Point> points);

        /* Returns the closest point in the tree to the query point, and its distance. */
        FindResult findBestMatch(const Point& query_point) const;

        constexpr const_iterator begin() const noexcept { return points_.begin(); }
        constexpr const_iterator end() const noexcept   { return points_.end(); }

        constexpr size_t size() const noexcept { return points_.size(); }
        constexpr bool empty() const noexcept { return points_.empty(); }

        constexpr const Matrix<double>& data() const noexcept { return points_; }

    private:
        Matrix<double> points_;
        std::vector<Node> nodes_;

        inline static constexpr size_t MAX_LEAF_ELEMENTS = 8;  /* The maximum number of points in a leaf node. */

        void buildTree();

        constexpr iterator node_begin(const Node& node) noexcept { return points_.begin() + node.first; }
        constexpr iterator node_end(const Node& node) noexcept { return points_.begin() + node.last; }

        constexpr const_iterator node_begin(const Node& node) const noexcept { return points_.begin() + node.first; }
        constexpr const_iterator node_end(const Node& node) const noexcept { return points_.begin() + node.last; }

        constexpr const_iterator node_cbegin(const Node& node) const noexcept { return points_.cbegin() + node.first; }
        constexpr const_iterator node_cend(const Node& node) const noexcept { return points_.cbegin() + node.last; }

        const Node& left_child(const Node& node) const noexcept { return nodes_[node.left]; }
        const Node& right_child(const Node& node) const noexcept { return nodes_[node.right]; }
    };

} // namespace gapp::detail

#endif // !GA_UTILITY_CONE_TREE_HPP