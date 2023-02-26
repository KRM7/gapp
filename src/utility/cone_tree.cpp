/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "cone_tree.hpp"
#include "functional.hpp"
#include "math.hpp"
#include "matrix.hpp"
#include <algorithm>
#include <numeric>
#include <functional>
#include <iterator>
#include <utility>
#include <cmath>
#include <cassert>

namespace genetic_algorithm::detail
{
    using iterator = ConeTree::iterator;
    using const_iterator = ConeTree::const_iterator;

    using Node = ConeTree::Node;
    using FindResult = ConeTree::FindResult;


    ConeTree::ConeTree(std::span<const Point> points)
    {
        if (points.empty()) return;

        points_ = Matrix(0, points[0].size(), 0.0);
        points_.reserve(points.size(), points[0].size());
        for (const Point& point : points) points_.append_row(point);

        nodes_.reserve(4 * points_.size() / MAX_LEAF_ELEMENTS);
        Node root{ .first = 0, .last = points_.size() };
        nodes_.push_back(root);

        buildTree();
    }

    /* Find the point in the range [first, last) furthest from a point (using Euclidean distances). */
    static inline const_iterator findFurthestElement(const_iterator first, const_iterator last, const_iterator from)
    {
        assert(std::distance(first, last) > 0);

        const_iterator furthest;
        double max_distance = -math::inf<double>;

        for (; first != last; ++first)
        {
            const double distance = math::euclideanDistanceSq(first->begin(), first->end(), from->begin());

            if (distance > max_distance)
            {
                furthest = first;
                max_distance = distance;
            }
        }

        return furthest;
    }

    /* Find the 2 partition points that will be used to split the range [first, last) into 2 parts. */
    static inline std::pair<const_iterator, const_iterator> partitionPoints(const_iterator first, const_iterator last)
    {
        assert(std::distance(first, last) > 0);

        const_iterator rand_elem = first;

        const_iterator first_elem = findFurthestElement(first, last, rand_elem);
        const_iterator second_elem = findFurthestElement(first, last, first_elem);

        return { first_elem, second_elem };
    }

    /* The center of the range of points [first, last) is the mean of the coords along each axis. */
    static inline Point findCenter(const_iterator first, const_iterator last)
    {
        assert(std::distance(first, last) > 0);

        const ptrdiff_t range_len = std::distance(first, last);

        Point center(first->begin(), first->end());

        for (++first; first != last; ++first)
        {
            std::transform(center.begin(), center.end(), first->begin(), center.begin(), std::plus{});
        }

        std::transform(center.begin(), center.end(), center.begin(), detail::divide_by(range_len));

        return center;
    }

    /* Find the Euclidean distance between the center point and the point in the range [first, last) furthest from it. */
    static inline double findRadius(const_iterator first, const_iterator last, const Point& center)
    {
        assert(std::distance(first, last) > 0);

        double max_distance = -math::inf<double>;

        for (; first != last; ++first)
        {
            const double distance = math::euclideanDistanceSq(center.begin(), center.end(), first->begin());
            max_distance = std::max(max_distance, distance);
        }

        return std::sqrt(max_distance);
    }

    /* Returns true if the node is a leaf node. */
    static inline bool isLeafNode(const Node& node) noexcept
    {
        return !static_cast<bool>(node.left & node.right);
    }

    /* Return the max possible inner product between the point and a point inside the node. */
    static inline double innerProductUpperBound(const Node& node, const Point& point, double point_norm)
    {
        const double center_prod = std::inner_product(point.begin(), point.end(), node.center.begin(), 0.0);

        return center_prod + point_norm * node.radius;
    }

    /* Find the best match in the range [first, last) using linear search. */
    static FindResult findBestMatchLinear(const Point& query_point, const_iterator first, const_iterator last)
    {
        assert(std::distance(first, last) > 0);
        assert(query_point.size() == first->size());

        FindResult best = { {}, -math::inf<double> };

        for (; first != last; ++first)
        {
            const double inner_prod = std::inner_product(query_point.begin(), query_point.end(), first->begin(), 0.0);

            if (inner_prod > best.prod)
            {
                best.elem = first;
                best.prod = inner_prod;
            }
        }

        return best;
    }

    void ConeTree::buildTree()
    {
        assert(nodes_.size() == 1);

        for (size_t i = 0; i < nodes_.size(); i++)
        {
            Node& node = nodes_[i];

            node.center = findCenter(node_cbegin(node), node_cend(node));
            node.radius = findRadius(node_cbegin(node), node_cend(node), node.center);

            /* Leaf node. */
            if (size_t(node.last - node.first) <= MAX_LEAF_ELEMENTS)
            {
                node.left = 0;
                node.right = 0;
            }
            /* Non-leaf node. */
            else
            {
                const auto partition_points = partitionPoints(node_cbegin(node), node_cend(node));

                auto middle = std::partition(node_begin(node), node_end(node), [&](auto point)
                {
                    const double left_dist = math::euclideanDistanceSq(partition_points.first->begin(), partition_points.first->end(), point.begin());
                    const double right_dist = math::euclideanDistanceSq(partition_points.second->begin(), partition_points.second->end(), point.begin());

                    return left_dist < right_dist;
                });

                /* Handle edge case where all of the points in [first, last) are the same (making sure both child ranges will be non-empty). */
                if (middle == node_begin(node)) ++middle;

                const size_t middle_idx = size_t(middle - points_.begin());

                Node left_child{ .first = node.first, .last = middle_idx };
                Node right_child{ .first = middle_idx, .last = node.last };

                nodes_.push_back(left_child);
                nodes_.push_back(right_child);

                node.left = nodes_.size() - 2;
                node.right = nodes_.size() - 1;
            }
        }
    }

    FindResult ConeTree::findBestMatch(const Point& query_point) const
    {
        assert(query_point.size() == points_.ncols());

        const double query_norm = math::euclideanNorm(query_point);

        static thread_local std::vector<const Node*> node_stack(nodes_.size() / 2);
        node_stack.clear();
        node_stack.push_back(&nodes_[0]);

        FindResult best{ {}, -math::inf<double> };

        while (!node_stack.empty())
        {
            const Node* node = node_stack.back();
            node_stack.pop_back();

            /* Skip node if it can't be better. */
            if (best.prod >= innerProductUpperBound(*node, query_point, query_norm)) continue;

            if (isLeafNode(*node))
            {
                auto [elem, inner_prod] = findBestMatchLinear(query_point, node_begin(*node), node_end(*node));
                if (inner_prod > best.prod)
                {
                    best.elem = elem;
                    best.prod = inner_prod;
                }
            }
            else
            {
                if (innerProductUpperBound(left_child(*node), query_point, query_norm) <
                    innerProductUpperBound(right_child(*node), query_point, query_norm))
                {
                    /* Visit right child first. */
                    node_stack.push_back(&left_child(*node));
                    node_stack.push_back(&right_child(*node));
                }
                else
                {
                    /* Visit left child first. */
                    node_stack.push_back(&right_child(*node));
                    node_stack.push_back(&left_child(*node));
                }
            }
        }

        return best;
    }

} // namespace genetic_algorithm::detail