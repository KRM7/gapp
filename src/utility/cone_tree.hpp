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
        struct FindResult
        {
            T* elem;
            double prod;
        };

        ConeTree() requires std::is_invocable_r_v<Point, Proj, T> = default;

        template<std::input_iterator Iter>
        requires std::is_invocable_r_v<Point, Proj, T>
        ConeTree(Iter first, Iter last, Proj proj = {});

        size_t size() const noexcept { return data_.size(); }
        
        std::vector<T>& data() noexcept { return data_; }
        const std::vector<T>& data() const noexcept { return data_; }

        /* Returns the closest point in the tree to the point, and its distance. */
        FindResult findBestMatch(const Point& query_point) const;

    private:
        struct Node
        {
            Point center;
            double radius;
            std::vector<T*> elements;

            std::unique_ptr<Node> left;
            std::unique_ptr<Node> right;
        };

        Node root_;
        Proj proj_;
        std::vector<T> data_;
        size_t dim_;

        /* The maximum number of elements allowed in a leaf node. */
        inline static constexpr size_t MAX_LEAF_ELEMENTS = 24;

        /* Find the element in a node furthest from an element (using Euclidean distances). */
        const T* findFurthestElement(const Node& node, const T& from) const;

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
        double innerProductUpperBound(const Point& point, double point_norm, const Node& node) const;

        /* Find the best match among node_elements using linear search. */
        FindResult findBestMatchLinear(const Point& query_point, const Node& node) const;

        /* Find the best match under the node. */
        FindResult findBestMatchDFS(const Point& query_point, const Node& node) const;
    };

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include <algorithm>
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
        data_(first, last),
        proj_(std::move(proj)),
        dim_(std::invoke(proj_, *first).size())
    {
        assert(std::all_of(first, last, [&](const auto& elem) { return std::invoke(proj_, elem).size() == dim_; }));

        root_.elements.reserve(data_.size());
        std::transform(data_.begin(), data_.end(), std::back_inserter(root_.elements), [](T& elem) { return &elem; });

        recursiveExpandNode(root_);
    }

    template<typename T, typename P>
    const T* ConeTree<T, P>::findFurthestElement(const Node& node, const T& from) const
    {
        assert(!node.elements.empty());

        const Point& point = std::invoke(proj_, from);

        const T* furthest = nullptr;
        double max_distance = -math::inf<double>;

        for (const T* elem : node.elements)
        {
            const Point& elem_point = std::invoke(proj_, *elem);
            const double distance = math::euclideanDistanceSq(elem_point, point);

            if (distance > max_distance)
            {
                furthest = elem;
                max_distance = distance;
            }
        }

        return furthest;
    }

    template<typename T, typename P>
    std::pair<const T*, const T*> ConeTree<T, P>::pickNodeSplitPoints(const Node& node) const
    {
        assert(!node.elements.empty());

        const T* rand_elem = node.elements.front(); 

        const T* first_elem  = findFurthestElement(node, *rand_elem);
        const T* second_elem = findFurthestElement(node, *first_elem);

        return { first_elem, second_elem };
    }

    template<typename T, typename P>
    std::vector<double> ConeTree<T, P>::findNodeCenter(const Node& node) const
    {
        assert(!node.elements.empty());

        const size_t ndim = std::invoke(proj_, node.elements.front()).size();

        Point center(ndim, 0.0);

        for (const T* elem : node.elements)
        {
            const Point& point = std::invoke(proj_, *elem);
            std::transform(center.begin(), center.end(), point.begin(), center.begin(), std::plus{});
        }

        std::transform(center.begin(), center.end(), center.begin(), detail::divide_by(double(node.elements.size())));

        return center;
    }

    template<typename T, typename P>
    double ConeTree<T, P>::findNodeRadius(const Node& node) const
    {
        assert(!node.elements.empty());

        double max_distance = -math::inf<double>;

        for (const T* elem : node.elements)
        {
            const Point& point = std::invoke(proj_, *elem);
            const double distance = math::euclideanDistanceSq(node.center, point);

            max_distance = std::max(max_distance, distance);
        }

        return std::sqrt(max_distance);
    }

    template<typename T, typename P>
    void ConeTree<T, P>::recursiveExpandNode(Node& node) const
    {
        assert(!node.elements.empty());

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

            const size_t child_size = std::max(1_sz, size_t(0.6 * node.elements.size()));
            node.left->elements.reserve(child_size);
            node.right->elements.reserve(child_size);

            const auto [left_elem, right_elem] = pickNodeSplitPoints(node);

            const Point& left_point = std::invoke(proj_, *left_elem);
            const Point& right_point = std::invoke(proj_, *right_elem);

            for (T* elem : node.elements)
            {
                const Point& elem_point = std::invoke(proj_, *elem);
                const double left_distance = math::euclideanDistanceSq(left_point, elem_point);
                const double right_distance = math::euclideanDistanceSq(right_point, elem_point);

                (left_distance <= right_distance) ?
                    node.left->elements.push_back(elem) :
                    node.right->elements.push_back(elem);
            }

            recursiveExpandNode(*node.left);
            recursiveExpandNode(*node.right);
        }
    }


    template<typename T, typename P>
    inline bool ConeTree<T, P>::isLeafNode(const Node& node) const noexcept
    {
        return (node.left == nullptr) && (node.right == nullptr);
    }

    template<typename T, typename P>
    inline double ConeTree<T, P>::innerProductUpperBound(const Point& point, double point_norm, const Node& node) const
    {
        const double center_prod = std::inner_product(point.begin(), point.end(), node.center.begin(), 0.0);
        
        return center_prod + point_norm * node.radius;
    }

    template<typename T, typename P>
    auto ConeTree<T, P>::findBestMatchLinear(const Point& query_point, const Node& node) const -> FindResult
    {
        assert(!node.elements.empty());
        assert(query_point.size() == dim_);

        FindResult best = { nullptr, -math::inf<double> };
        
        for (T* elem : node.elements)
        {
            const Point& node_point = std::invoke(proj_, *elem);
            const double inner_prod = std::inner_product(query_point.begin(), query_point.end(), node_point.begin(), 0.0);

            if (inner_prod > best.prod)
            {
                best.elem = elem;
                best.prod = inner_prod;
            }
        }
        
        return best;
    }

    template<typename T, typename P>
    auto ConeTree<T, P>::findBestMatchDFS(const Point& query_point, const Node& root) const -> FindResult
    {
        assert(!root.elements.empty());
        assert(query_point.size() == dim_);

        const double query_norm = math::euclideanNorm(query_point);

        std::vector<const Node*> node_stack;
        node_stack.reserve(size());
        node_stack.push_back(&root);

        FindResult best{ nullptr, -math::inf<double> };

        while (!node_stack.empty())
        {
            const Node* cur_node = node_stack.back();
            node_stack.pop_back();

            /* Skip node if it can't be better. */
            if (best.prod >= innerProductUpperBound(query_point, query_norm, *cur_node)) continue;
             
            if (isLeafNode(*cur_node))
            {
                auto [elem, inner_prod] = findBestMatchLinear(query_point, *cur_node);
                if (inner_prod > best.prod)
                {
                    best.elem = elem;
                    best.prod = inner_prod;
                }
            }
            else if (cur_node->left == nullptr)
            {
                /* Only need to visit right child. */
                node_stack.push_back(cur_node->right.get());
            }
            else if (cur_node->right == nullptr)
            {
                /* Only need to visit left child. */
                node_stack.push_back(cur_node->left.get());
            }
            else
            {
                if (innerProductUpperBound(query_point, query_norm, *(cur_node->left)) <
                    innerProductUpperBound(query_point, query_norm, *(cur_node->right)))
                {
                    /* Visit right child first. */
                    node_stack.push_back(cur_node->left.get());
                    node_stack.push_back(cur_node->right.get());
                }
                else
                {
                    /* Visit left child first. */
                    node_stack.push_back(cur_node->right.get());
                    node_stack.push_back(cur_node->left.get());
                }
            }
        }

        return best;
    }

    template<typename T, typename P>
    inline auto ConeTree<T, P>::findBestMatch(const Point& point) const -> FindResult
    {
        return findBestMatchDFS(point, root_);
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_CONE_TREE_HPP