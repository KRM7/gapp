/* Example showing the usage of metrics in the GAs. */

#include "gapp.hpp"
#include <vector>
#include <iostream>
#include <format>
#include <cassert>

using namespace gapp;

struct MyMetric : public metrics::Monitor<MyMetric, std::vector<double>>
{
    // track the fitness of the first solution in the population
    void update(const GaInfo& ga) override { data_.push_back(ga.fitness_matrix()[0][0]); }
};

int main()
{
    RCGA GA;

    // set the metrics to track and run

    GA.track(metrics::FitnessMin{}, metrics::FitnessMax{}, MyMetric{});
    GA.solve(problems::Sphere{ 10 }, Bounds{ -5.0, 5.0 });

    // accessing the recorded metric values

    const MyMetric& metric = GA.get_metric<MyMetric>();

    std::cout << "The values of MyMetric throughout the run:\n";
    for (size_t gen = 0; gen < metric.size(); gen++)
    {
        std::cout << std::format("Generation {}\t| {:.6f}\n", gen + 1, metric[gen]);
    }

     // trying to read an untracked metric

    [[maybe_unused]] const auto* hypervol = GA.get_metric_if<metrics::AutoHypervolume>();
    assert(hypervol == nullptr);
}
