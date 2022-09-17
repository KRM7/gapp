/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_CONE_TREE_HPP
#define GA_UTILITY_CONE_TREE_HPP

#include <vector>
#include <functional>
#include <concepts>
#include <cstddef>
#include <utility>
#include "math.hpp"

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
    template<typename T, typename Proj = std::identity>
    class ConeTree
    {
    public:
        using iterator = std::vector<T>::iterator;
        using const_iterator = std::vector<T>::const_iterator;

        struct FindResult
        {
            iterator elem;
            double prod;
        };

        ConeTree() requires std::is_invocable_r_v<Point, Proj, T> = default;

        template<std::input_iterator Iter>
        requires std::is_invocable_r_v<Point, Proj, T>
        ConeTree(Iter first, Iter last, Proj proj = {});


        /* Returns the closest point in the tree to the point, and its distance. */
        FindResult findBestMatch(const Point& query_point) const;


        size_t size() const noexcept { return points_.size(); }

        std::vector<T>& data() noexcept             { return points_; }
        const std::vector<T>& data() const noexcept { return points_; }

        iterator begin()                { return points_.begin(); }
        iterator end()                  { return points_.end(); }
        const_iterator begin() const    { return points_.begin(); }
        const_iterator end() const      { return points_.end(); }

    private:
        struct Node
        {
            Point center;
            double radius;
            iterator first;
            iterator last;

            size_t left;
            size_t right;
        };

        std::vector<T> points_;
        std::vector<Node> nodes_;
        Proj proj_;
        size_t dim_;

        /* The maximum number of elements allowed in a leaf node. */
        inline static constexpr size_t MAX_LEAF_ELEMENTS = 22;


        /* Find the element in the range [first, last) furthest from an element (using Euclidean distances). */
        const_iterator findFurthestElement(const_iterator first, const_iterator last, const_iterator from) const;

        /* Find the 2 partition points that will be used to split the range [first, last) into 2 parts. */
        std::pair<const_iterator, const_iterator> partitionPoints(const_iterator first, const_iterator last) const;

        /* The center of the range [first, last) is the mean of the elements along each axis. */
        Point findCenter(const_iterator first, const_iterator last) const;

        /* Find the Euclidean distance between the center and the element in the range [first, last) furthest from it. */
        double findRadius(const_iterator first, const_iterator last, const Point& center) const;

        /* Expand the tree from the root. */
        void buildTree();


        /* Returns true if the node is a leaf node. */
        bool isLeafNode(const Node& node) const noexcept;

        /* Return the max possible inner product between the point and a point inside the node. */
        double innerProductUpperBound(const Point& point, double point_norm, const Node& node) const;

        /* Find the best match in the range [first, last) using linear search. */
        FindResult findBestMatchLinear(const Point& query_point, iterator first, iterator last) const;

        /* Find the best match under the node. */
        FindResult findBestMatchDFS(const Point& query_point, const Node& node) const;
    };

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include <algorithm>
#include <iterator>
#include <numeric>
#include <limits>
#include "algorithm.hpp"
#include "functional.hpp"
#include "utility.hpp"

namespace genetic_algorithm::detail
{
    template<typename Iter, typename Proj>
    ConeTree(Iter, Iter, Proj) -> ConeTree<typename std::iterator_traits<Iter>::value_type, Proj>;

    template<typename T, typename P>
    template<std::input_iterator Iter>
    requires std::is_invocable_r_v<Point, P, T>
    ConeTree<T, P>::ConeTree(Iter first, Iter last, P proj) :
        points_(first, last),
        proj_(std::move(proj))
    {
        dim_ = std::invoke(proj_, points_.front()).size();
        nodes_.reserve(4 * points_.size() / MAX_LEAF_ELEMENTS);

        Node root{ .first = points_.begin(), .last = points_.end() };
        nodes_.push_back(root);

        buildTree();
    }

    template<typename T, typename P>
    auto ConeTree<T, P>::findFurthestElement(const_iterator first, const_iterator last, const_iterator from) const -> const_iterator
    {
        assert(std::distance(first, last) > 0);

        const Point& from_point = std::invoke(proj_, *from);

        const_iterator furthest;
        double max_distance = -math::inf<double>;

        for (; first != last; ++first)
        {
            const Point& this_point = std::invoke(proj_, *first);
            const double distance = math::euclideanDistanceSq(this_point, from_point);

            if (distance > max_distance)
            {
                furthest = first;
                max_distance = distance;
            }
        }

        return furthest;
    }

    template<typename T, typename P>
    inline auto ConeTree<T, P>::partitionPoints(const_iterator first, const_iterator last) const -> std::pair<const_iterator, const_iterator>
    {
        assert(std::distance(first, last) > 0);

        const_iterator rand = first;

        const_iterator first_elem = findFurthestElement(first, last, rand);
        const_iterator second_elem = findFurthestElement(first, last, first_elem);

        return { first_elem, second_elem };
    }

    template<typename T, typename P>
    Point ConeTree<T, P>::findCenter(const_iterator first, const_iterator last) const
    {
        assert(std::distance(first, last) > 0);

        const size_t ndim = std::invoke(proj_, *first).size();
        const ptrdiff_t nrange = std::distance(first, last);

        Point center(ndim, 0.0);

        for (; first != last; ++first)
        {
            const Point& point = std::invoke(proj_, *first);
            std::transform(center.begin(), center.end(), point.begin(), center.begin(), std::plus{});
        }

        std::transform(center.begin(), center.end(), center.begin(), detail::divide_by(nrange));

        return center;
    }

    template<typename T, typename P>
    double ConeTree<T, P>::findRadius(const_iterator first, const_iterator last, const Point& center) const
    {
        assert(std::distance(first, last) > 0);

        double max_distance = -math::inf<double>;

        for (; first != last; ++first)
        {
            const Point& point = std::invoke(proj_, *first);
            const double distance = math::euclideanDistanceSq(center, point);

            max_distance = std::max(max_distance, distance);
        }

        return std::sqrt(max_distance);
    }

    template<typename T, typename P>
    void ConeTree<T, P>::buildTree()
    {
        assert(nodes_.size() == 1);

        for (size_t i = 0; i < nodes_.size(); i++)
        {
            Node& node = nodes_[i];

            node.center = findCenter(node.first, node.last);
            node.radius = findRadius(node.first, node.last, node.center);

            /* Leaf node. */
            if (size_t(node.last - node.first) <= MAX_LEAF_ELEMENTS)
            {
                node.left = 0;
                node.right = 0;
            }
            /* Non-leaf node. */
            else
            {
                const auto [left_elem, right_elem] = partitionPoints(node.first, node.last);

                const Point left_point = std::invoke(proj_, *left_elem);
                const Point right_point = std::invoke(proj_, *right_elem);

                auto middle = std::partition(node.first, node.last, [&](const T& elem)
                {
                    const Point& this_point = std::invoke(proj_, elem);

                    const double left_dist = math::euclideanDistanceSq(left_point, this_point);
                    const double right_dist = math::euclideanDistanceSq(right_point, this_point);

                    return left_dist < right_dist;
                });

                /* Handle edge case where all of the points in [first, last) are the same (make sure both child ranges will be non-empty). */
                if (middle == node.first) ++middle;

                Node left{ .first = node.first, .last = middle };
                Node right{ .first = middle, .last = node.last };

                node.left = nodes_.size();
                node.right = nodes_.size() + 1;

                nodes_.push_back(left);
                nodes_.push_back(right);
            }
        }
    }


    template<typename T, typename P>
    inline bool ConeTree<T, P>::isLeafNode(const Node& node) const noexcept
    {
        //return (node.left == nullptr) && (node.right == nullptr);
        return !(node.left & node.right);
    }

    template<typename T, typename P>
    inline double ConeTree<T, P>::innerProductUpperBound(const Point& point, double point_norm, const Node& node) const
    {
        const double center_prod = std::inner_product(point.begin(), point.end(), node.center.begin(), 0.0);
        
        return center_prod + point_norm * node.radius;
    }

    template<typename T, typename P>
    auto ConeTree<T, P>::findBestMatchLinear(const Point& query_point, iterator first, iterator last) const -> FindResult
    {
        assert(std::distance(first, last) > 0);
        assert(query_point.size() == dim_);

        FindResult best = { {}, -math::inf<double> };
        
        for (; first != last; ++first)
        {
            const Point& this_point = std::invoke(proj_, *first);
            const double inner_prod = std::inner_product(query_point.begin(), query_point.end(), this_point.begin(), 0.0);

            if (inner_prod > best.prod)
            {
                best.elem = first;
                best.prod = inner_prod;
            }
        }
        
        return best;
    }

    template<typename T, typename P>
    auto ConeTree<T, P>::findBestMatchDFS(const Point& query_point, const Node& root) const -> FindResult
    {
        assert(std::distance(root.first, root.last) > 0);
        assert(query_point.size() == dim_);

        const double query_norm = math::euclideanNorm(query_point);

        static thread_local std::vector<const Node*> node_stack(nodes_.size() / 2);
        node_stack.clear();
        node_stack.push_back(&root);

        FindResult best{ {}, -math::inf<double> };

        while (!node_stack.empty())
        {
            const Node* cur_node = node_stack.back();
            node_stack.pop_back();

            /* Skip node if it can't be better. */
            if (best.prod >= innerProductUpperBound(query_point, query_norm, *cur_node)) continue;
             
            if (isLeafNode(*cur_node))
            {
                auto [elem, inner_prod] = findBestMatchLinear(query_point, cur_node->first, cur_node->last);
                if (inner_prod > best.prod)
                {
                    best.elem = elem;
                    best.prod = inner_prod;
                }
            }
            else
            {
                if (innerProductUpperBound(query_point, query_norm, nodes_[cur_node->left]) <
                    innerProductUpperBound(query_point, query_norm, nodes_[cur_node->right]))
                {
                    /* Visit right child first. */
                    node_stack.push_back(&nodes_[cur_node->left]);
                    node_stack.push_back(&nodes_[cur_node->right]);
                }
                else
                {
                    /* Visit left child first. */
                    node_stack.push_back(&nodes_[cur_node->right]);
                    node_stack.push_back(&nodes_[cur_node->left]);
                }
            }
        }

        return best;
    }

    template<typename T, typename P>
    inline auto ConeTree<T, P>::findBestMatch(const Point& point) const -> FindResult
    {
        return findBestMatchDFS(point, nodes_[0]);
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_CONE_TREE_HPP