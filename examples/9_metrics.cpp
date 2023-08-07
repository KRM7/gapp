/* Example showing the usage of metrics in the GAs. */

#include "gapp.hpp"
#include <vector>
#include <iostream>
#include <format>
#include <cassert>

using namespace gapp;

struct MyMetric : public metrics::Monitor<MyMetric, std::vector<double>>
{
    double value_at(size_t generation) const noexcept { return data_[generation]; }
    void initialize(const GaInfo&) override { data_.clear(); }
    void update(const GaInfo& ga) override { data_.push_back(ga.fitness_matrix()[0][0]); }
};

int main()
{
    RCGA GA;
    GA.track(metrics::FitnessMin{}, metrics::FitnessMax{}, MyMetric{});
    GA.solve(problems::Sphere{ 10 }, Bounds{ -5.0, 5.0 });

    const MyMetric& metric = GA.get_metric<MyMetric>();

    std::cout << "The values of MyMetric throughout the run:\n";
    for (size_t gen = 0; gen < metric.size(); gen++)
    {
        std::cout << std::format("Generation {}\t| {:.6f}\n", gen + 1, metric[gen]);
    }

    const auto* hypervol = GA.get_metric_if<metrics::AutoHypervolume>(); // untracked metric
    assert(hypervol == nullptr);
}
