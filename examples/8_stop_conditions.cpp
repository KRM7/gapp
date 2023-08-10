/* Example showing the usage of the stop conditions in the GAs. */

#include "gapp.hpp"
#include <vector>
#include <iostream>
#include <format>

using namespace gapp;

class MyStopCondition : public stopping::StopCondition
{
    void initialize(const GaInfo&) override {}

    bool stop_condition(const GaInfo& ga) override
    {
        return ga.num_fitness_evals() >= 4000;
    }
};

int main()
{
    BinaryGA GA;

    GA.solve(problems::Sphere{ 10 });
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);

    GA.solve(problems::Sphere{ 10 }, /* generations */ 375);
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);

    GA.max_gen(755);
    GA.solve(problems::Sphere{ 10 });
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);

    GA.solve(problems::Sphere{ 10 }, /* generations */ 175);
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);

    // early-stop conditions

    GA.stop_condition(stopping::FitnessBestStall{ 3 });
    GA.solve(problems::Sphere{ 10 });
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);

    GA.stop_condition(nullptr); // same as GA.stop_condition(stopping::NoEarlyStop{});
    GA.solve(problems::Sphere{ 10 });
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);

    // composite early-stop conditions

    GA.stop_condition(stopping::FitnessBestStall{ 2 } && stopping::FitnessMeanStall{ 3 });
    GA.solve(problems::Sphere{ 10 }, 5000);
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);

    // custom early-stop conditions

    GA.stop_condition([](const GaInfo& ga) { return ga.num_fitness_evals() >= 10000; });
    GA.solve(problems::Sphere{ 10 });
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);

    GA.stop_condition(MyStopCondition{});
    GA.solve(problems::Sphere{ 10 });
    std::cout << std::format("The GA ran for {} generations.\n", GA.generation_cntr() + 1);
}
