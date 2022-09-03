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

    private:
        using Point = std::vector<double>;

        struct Node
        {
            Point center;
            double radius_sq;
            std::vector<T*> elements;

            std::unique_ptr<Node> left;
            std::unique_ptr<Node> right;
        };

        Node root_;
        Proj proj_;
        std::vector<T> data_;

        inline static constexpr size_t MAX_LEAF_ELEMENTS = 1;   /* The maximum number of elements allowed in a leaf node. */

        /* Find the 2 center points that will be used to split the node elements into 2 parts. */
        std::pair<T*, T*> splitData(const std::vector<T*>& node_elements) const;

        /* The center of the node is the mean of the element points along each axis. */
        std::vector<double> findNodeCenter(const std::vector<T*>& node_elements) const;

        /* Find the square of the euclidean distance between the center of the node, and the element furthest from the center. */
        double findNodeRadius(const std::vector<T*>& node_elements, const Point& center) const;

        /* Expand the tree starting from a node. */
        void buildTree(Node& node) const;
    };

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include <algorithm>
#include <numeric>
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
        data_(first, last), proj_(std::forward<Projection>(proj))
    {
        assert(std::all_of(first, last, [&](const auto& elem) { return std::invoke(proj, elem).size() == std::invoke(proj, *first).size(); }));

        std::transform(data_.begin(), data_.end(), std::back_inserter(root_.elements), [](T& elem) { return &elem; });

        buildTree(root_);
    }

    template<typename T, typename P>
    std::pair<T*, T*> ConeTree<T, P>::splitData(const std::vector<T*>& node_elements) const
    {
        const T* rand_node = node_elements.front();

        const auto distances1 = detail::map(node_elements, [&](const T* elem)
        {
            const Point& p1 = std::invoke(proj_, *rand_node);
            const Point& p2 = std::invoke(proj_, *elem);

            return detail::euclideanDistanceSq(p1, p2);
        });

        const size_t first_idx = detail::argmax(distances1.begin(), distances1.begin(), distances1.end());

        const auto distances2 = detail::map(node_elements, [&](const T* elem)
        {
            const Point& p1 = std::invoke(proj_, *node_elements[first_idx]);
            const Point& p2 = std::invoke(proj_, *elem);

            return detail::euclideanDistanceSq(p1, p2);
        });

        const size_t second_idx = detail::argmax(distances2.begin(), distances2.begin(), distances2.end());

        return { node_elements[first_idx], node_elements[second_idx] };
    }

    template<typename T, typename P>
    std::vector<double> ConeTree<T, P>::findNodeCenter(const std::vector<T*>& node_elements) const
    {
        const size_t ndim = std::invoke(proj_, node_elements.front()).size();

        Point center = std::accumulate(node_elements.begin(), node_elements.end(), Point(ndim, 0.0),
        [this](Point acc, const T* elem)
        {
            const Point& point = std::invoke(proj_, *elem);
            std::transform(acc.begin(), acc.end(), point.begin(), acc.begin(), std::plus{});

            return std::move(acc);
        });

        std::transform(center.begin(), center.end(), center.begin(), detail::divide_by(double(node_elements.size())));

        return center;
    }

    template<typename T, typename P>
    double ConeTree<T, P>::findNodeRadius(const std::vector<T*>& node_elements, const Point& center) const
    {
        const auto distances = detail::map(node_elements, [&](const T* elem)
        {
            const Point& point = std::invoke(proj_, *elem);
            return detail::euclideanDistanceSq(center, point);
        });

        const size_t furthest_idx = detail::argmax(distances.begin(), distances.begin(), distances.end());
        const Point& furthest_point = std::invoke(proj_, *node_elements[furthest_idx]);

        return detail::euclideanDistanceSq(center, furthest_point);
    }

    template<typename T, typename P>
    void ConeTree<T, P>::buildTree(Node& node) const
    {
        node.center = findNodeCenter(node.elements);
        node.radius_sq = findNodeRadius(node.elements, node.center);

        /* Leaf node. */
        if (node.elements.size() <= MAX_LEAF_ELEMENTS)
        {
            node.left = nullptr;
            node.right = nullptr;
        }
        /* Non-leaf node. */
        else
        {
            const auto [left, right] = splitData(node.elements);

            const Point& left_point = std::invoke(proj_, *left);
            const Point& right_point = std::invoke(proj_, *right);

            const auto distances = detail::map(node.elements, [&](const T* elem)
            {
                const Point& point = std::invoke(proj_, *elem);
                const double left_distance = detail::euclideanDistanceSq(left_point, point);
                const double right_distance = detail::euclideanDistanceSq(right_point, point);

                return std::make_pair(left_distance, right_distance);
            });

            node.left = std::make_unique<Node>();
            node.right = std::make_unique<Node>();

            for (size_t i = 0; i < distances.size(); i++)
            {
                (distances[i].first <= distances[i].second) ?
                    node.left->elements.push_back(node.elements[i]) :
                    node.right->elements.push_back(node.elements[i]);
            }

            buildTree(*node.left);
            buildTree(*node.right);
        }
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_CONE_TREE_HPP