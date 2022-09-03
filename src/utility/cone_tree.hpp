/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_CONE_TREE_HPP
#define GA_UTILITY_CONE_TREE_HPP

#include <vector>
#include <functional>
#include <concepts>
#include <cstddef>
#include <memory>
#include <utility>
#include "math.hpp"

namespace genetic_algorithm::detail
{
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
        template<std::input_iterator Iter, typename P = std::identity>
        requires std::is_invocable_r_v<Point, P, T>
        ConeTree(Iter first, Iter last, P&& p = std::identity{});

        size_t size() const noexcept { return data_.size(); }

        /* Returns the closest point in the tree to the point, and its distance. */
        std::pair<const T*, double> findBestMatch(const Point& query_point) const;

    private:
        using Point = std::vector<double>;

        struct Node
        {
            Point center;
            double radius;
            std::vector<const T*> elements;

            std::unique_ptr<Node> left;
            std::unique_ptr<Node> right;
        };

        Node root_;
        Proj proj_;
        std::vector<T> data_;
        size_t dim_;

        /* The maximum number of elements allowed in a leaf node. */
        inline static constexpr size_t MAX_LEAF_ELEMENTS = 4;


        /* Find the 2 center points that will be used to split the node elements into 2 parts. */
        std::pair<const T*, const T*> pickNodeSplitPoints(const Node& node) const;

        /* The center of the node is the mean of the element points along each axis. */
        std::vector<double> findNodeCenter(const Node& node_elements) const;

        /* Find the square of the euclidean distance between the center of the node, and the element furthest from the center. */
        double findNodeRadius(const Node& node) const;

        /* Expand the tree starting from a node. */
        void recursiveExpandNode(Node& node) const;


        /* Returns true if the node is a leaf node. */
        bool isLeafNode(const Node& node) const noexcept;

        /* Return the max possible inner product between the point and a point inside the node. */
        double innerProductUpperBound(const Point& point, const Node& node) const;

        /* Find the best match among node_elements using linear search. */
        std::pair<const T*, double> findBestMatchLinear(const Point& query_point, const Node& node) const;

        /* Find the best match under the node. */
        std::pair<const T*, double> findBestMatchDFS(const Point& query_point, const Node& node) const;
    };

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include <stack>
#include <algorithm>
#include <numeric>
#include <limits>
#include "algorithm.hpp"
#include "functional.hpp"
#include "utility.hpp"

namespace genetic_algorithm::detail
{
    template<typename Iter, typename Proj>
    ConeTree(Iter, Iter, Proj&&) -> ConeTree<typename std::iterator_traits<Iter>::value_type, Proj>;


    template<typename T, typename P>
    template<std::input_iterator Iter, typename Projection>
    requires std::is_invocable_r_v<Point, Projection, T>
    ConeTree<T, P>::ConeTree(Iter first, Iter last, Projection&& proj) :
        data_(first, last),
        proj_(std::forward<Projection>(proj)),
        dim_(std::invoke(proj_, *first).size())
    {
        assert(std::all_of(first, last, [&](const auto& elem) { return std::invoke(proj_, elem).size() == dim_; }));

        root_.elements.reserve(data_.size());
        std::transform(data_.begin(), data_.end(), std::back_inserter(root_.elements), [](const T& elem) { return &elem; });

        recursiveExpandNode(root_);
    }

    template<typename T, typename P>
    std::pair<const T*, const T*> ConeTree<T, P>::pickNodeSplitPoints(const Node& node) const
    {
        assert(!node.elements.empty());

        const T* rand_node = node.elements.front();

        auto distances = detail::map(node.elements, [&](const T* elem)
        {
            const Point& p1 = std::invoke(proj_, *rand_node);
            const Point& p2 = std::invoke(proj_, *elem);

            return detail::euclideanDistanceSq(p1, p2);
        });

        const size_t first_idx = detail::argmax(distances.begin(), distances.begin(), distances.end());

        std::transform(node.elements.begin(), node.elements.end(), distances.begin(),
        [&](const T* elem)
        {
            const Point& p1 = std::invoke(proj_, *node.elements[first_idx]);
            const Point& p2 = std::invoke(proj_, *elem);

            return detail::euclideanDistanceSq(p1, p2);
        });

        const size_t second_idx = detail::argmax(distances.begin(), distances.begin(), distances.end());

        return { node.elements[first_idx], node.elements[second_idx] };
    }

    template<typename T, typename P>
    std::vector<double> ConeTree<T, P>::findNodeCenter(const Node& node) const
    {
        assert(!node.elements.empty());

        const size_t ndim = std::invoke(proj_, node.elements.front()).size();

        Point center = std::accumulate(node.elements.begin(), node.elements.end(), Point(ndim, 0.0),
        [&](Point acc, const T* elem)
        {
            const Point& point = std::invoke(proj_, *elem);
            std::transform(acc.begin(), acc.end(), point.begin(), acc.begin(), std::plus{});

            return std::move(acc);
        });

        std::transform(center.begin(), center.end(), center.begin(), detail::divide_by(double(node.elements.size())));

        return center;
    }

    template<typename T, typename P>
    double ConeTree<T, P>::findNodeRadius(const Node& node) const
    {
        const auto distances = detail::map(node.elements, [&](const T* elem)
        {
            const Point& point = std::invoke(proj_, *elem);
            return detail::euclideanDistanceSq(node.center, point);
        });

        const size_t furthest_idx = detail::argmax(distances.begin(), distances.begin(), distances.end());

        return std::sqrt(distances[furthest_idx]);
    }

    template<typename T, typename P>
    void ConeTree<T, P>::recursiveExpandNode(Node& node) const
    {
        node.center = findNodeCenter(node);
        node.radius = findNodeRadius(node);

        /* Leaf node. */
        if (node.elements.size() <= MAX_LEAF_ELEMENTS)
        {
            node.left = nullptr;
            node.right = nullptr;
        }
        /* Non-leaf node. */
        else
        {
            node.left = std::make_unique<Node>();
            node.right = std::make_unique<Node>();

            node.left->elements.reserve(size_t(0.6 * node.elements.size()));
            node.right->elements.reserve(size_t(0.6 * node.elements.size()));

            const auto [left_elem, right_elem] = pickNodeSplitPoints(node);

            const Point& left_point = std::invoke(proj_, *left_elem);
            const Point& right_point = std::invoke(proj_, *right_elem);

            const auto distances = detail::map(node.elements, [&](const T* elem)
            {
                const Point& elem_point = std::invoke(proj_, *elem);
                const double left_distance = detail::euclideanDistanceSq(left_point, elem_point);
                const double right_distance = detail::euclideanDistanceSq(right_point, elem_point);

                return std::make_pair(left_distance, right_distance);
            });

            for (size_t idx = 0; idx < distances.size(); idx++)
            {
                (distances[idx].first <= distances[idx].second) ?
                    node.left->elements.push_back(node.elements[idx]) :
                    node.right->elements.push_back(node.elements[idx]);
            }

            recursiveExpandNode(*node.left);
            recursiveExpandNode(*node.right);
        }
    }


    template<typename T, typename P>
    bool ConeTree<T, P>::isLeafNode(const Node& node) const noexcept
    {
        return (node.left == nullptr) && (node.right == nullptr);
    }

    template<typename T, typename P>
    double ConeTree<T, P>::innerProductUpperBound(const Point& point, const Node& node) const
    {
        const double center_prod = std::inner_product(point.begin(), point.end(), node.center.begin(), 0.0);
        
        return center_prod + detail::euclideanNorm(point) * node.radius;
    }

    template<typename T, typename P>
    std::pair<const T*, double> ConeTree<T, P>::findBestMatchLinear(const Point& query_point, const Node& node) const
    {
        assert(!node.elements.empty());
        assert(query_point.size() == dim_);

        const auto products = detail::map(node.elements, [&](const T* elem)
        {
            const Point& node_point = std::invoke(proj_, *elem);
            return std::inner_product(query_point.begin(), query_point.end(), node_point.begin(), 0.0);
        });

        const size_t best_idx = detail::argmax(products.begin(), products.begin(), products.end());

        return { node.elements[best_idx], products[best_idx] };
    }

    template<typename T, typename P>
    std::pair<const T*, double> ConeTree<T, P>::findBestMatchDFS(const Point& query_point, const Node& root) const
    {
        assert(!root.elements.empty());
        assert(query_point.size() == dim_);

        std::stack<const Node*> nodes;
        nodes.push(&root);

        const T* best_match = nullptr;
        double max_prod = -std::numeric_limits<double>::infinity();

        while (!nodes.empty())
        {
            const Node* cur_node = nodes.top();
            nodes.pop();

            /* Skip node if it can't be better. */
            if (max_prod >= innerProductUpperBound(query_point, *cur_node)) continue;
             
            if (isLeafNode(*cur_node))
            {
                auto [elem, inner_prod] = findBestMatchLinear(query_point, *cur_node);
                if (inner_prod > max_prod)
                {
                    best_match = elem;
                    max_prod = inner_prod;
                }
            }
            else if (cur_node->left == nullptr)
            {
                /* Only need to visit right child. */
                nodes.push(cur_node->right.get());
            }
            else if (cur_node->right == nullptr)
            {
                /* Only need to visit left child. */
                nodes.push(cur_node->left.get());
            }
            else
            {
                if (innerProductUpperBound(query_point, *(cur_node->left)) < innerProductUpperBound(query_point, *(cur_node->right)))
                {
                    /* Visit right child first. */
                    nodes.push(cur_node->left.get());
                    nodes.push(cur_node->right.get());
                }
                else
                {
                    /* Visit left child first. */
                    nodes.push(cur_node->right.get());
                    nodes.push(cur_node->left.get());
                }
            }
        }

        return { best_match, max_prod };
    }

    template<typename T, typename P>
    std::pair<const T*, double> ConeTree<T, P>::findBestMatch(const Point& point) const
    {
        return findBestMatchDFS(point, root_);
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_CONE_TREE_HPP