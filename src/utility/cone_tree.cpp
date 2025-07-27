/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "cone_tree.hpp"
#include "algorithm.hpp"
#include "functional.hpp"
#include "math.hpp"
#include "matrix.hpp"
#include "utility.hpp"
#include <algorithm>
#include <numeric>
#include <functional>
#include <iterator>
#include <utility>
#include <cmath>

namespace gapp::detail
{
    using iterator       = ConeTree::iterator;
    using const_iterator = ConeTree::const_iterator;

    using Point    = ConeTree::Point;
    using PointRef = ConeTree::PointRef;

    using Node       = ConeTree::Node;
    using FindResult = ConeTree::FindResult;


    ConeTree::ConeTree(std::span<const Point> points)
    {
        if (points.empty()) return;

        points_ = Matrix(0, points[0].size(), 0.0);
        points_.reserve(points.size(), points[0].size());
        for (const Point& point : points) points_.append_row(point);

        nodes_.reserve(4 * points_.size() / MAX_LEAF_ELEMENTS);
        nodes_.push_back({ .first = 0, .last = points_.size() });

        buildTree();
    }

    /* Find the point in the range [first, last) furthest from a point (using Euclidean distances). */
    static inline const_iterator findFurthestElement(const_iterator first, const_iterator last, const_iterator from)
    {
        return detail::max_element(first, last, std::bind_front(math::euclideanDistanceSq, *from));
    }

    /* Find the 2 partition points that will be used to split the range [first, last) into 2 parts. */
    static inline std::pair<const_iterator, const_iterator> partitionPoints(const_iterator first, const_iterator last)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);

        const_iterator rand_elem = first;

        const_iterator first_elem = findFurthestElement(first, last, rand_elem);
        const_iterator second_elem = findFurthestElement(first, last, first_elem);

        return { first_elem, second_elem };
    }

    /* The center of the range of points [first, last) is the mean of the coords along each axis. */
    static inline Point findCenter(const_iterator first, const_iterator last)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);

        const ptrdiff_t range_len = std::distance(first, last);

        Point center(*first);

        for (++first; first != last; ++first)
        {
            std::transform(center.begin(), center.end(), first->begin(), center.begin(), std::plus{});
        }
        std::transform(center.begin(), center.end(), center.begin(), detail::divide_by(range_len));

        return center;
    }

    /* Find the Euclidean distance between the center point and the point in the range [first, last) furthest from it. */
    static inline double findRadius(const_iterator first, const_iterator last, PointRef center)
    {
        auto distance = std::bind_front(math::euclideanDistanceSq, center);
        auto furthest = detail::max_element(first, last, distance);

        return std::sqrt(distance(*furthest));
    }

    /* Returns true if the node is a leaf node. */
    static inline bool isLeafNode(const Node& node) noexcept
    {
        return (node.left == 0) && (node.right == 0);
    }

    /* Return the max possible inner product between the point and a point inside the node. */
    static inline double innerProductUpperBound(const Node& node, PointRef point, double point_norm)
    {
        const double center_prod = std::inner_product(point.begin(), point.end(), node.center.begin(), 0.0);

        return center_prod + point_norm * node.radius;
    }

    /* Find the best match in the range [first, last) using linear search. */
    static FindResult findBestMatchLinear(PointRef query_point, const_iterator first, const_iterator last)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);
        GAPP_ASSERT(query_point.size() == first->size());

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
        GAPP_ASSERT(nodes_.size() == 1);

        for (size_t i = 0; i < nodes_.size(); i++)
        {
            Node& node = nodes_[i];

            node.center = findCenter(node_cbegin(node), node_cend(node));
            node.radius = findRadius(node_cbegin(node), node_cend(node), node.center);

            /* Leaf node. */
            if (node.last - node.first <= MAX_LEAF_ELEMENTS)
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
                    const double left_dist = math::euclideanDistanceSq(point, *partition_points.first);
                    const double right_dist = math::euclideanDistanceSq(point, *partition_points.second);

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

    FindResult ConeTree::findBestMatch(PointRef query_point) const
    {
        if (points_.empty()) return { points_.end(), 0.0 };

        GAPP_ASSERT(query_point.size() == points_.ncols());

        const double query_norm = math::euclideanNorm(query_point);

        thread_local std::vector<const Node*> node_stack(nodes_.size() / 2);
        node_stack.clear();
        node_stack.push_back(&nodes_.front());

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
                GAPP_ASSERT(node->left && node->right);

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

} // namespace gapp::detail
