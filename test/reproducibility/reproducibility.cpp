/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include "gapp.hpp"
#include <iostream>

using namespace gapp;

int main()
{
    rng::prng.seed(0x3da99432ab975d26LL);
    std::cout << "integer-1: " << rng::randomInt(0, 100) << "\n";

    rng::prng.seed(0x3da99432ab975d26LL);
    std::cout << "integer-2: " << rng::randomInt(0, 100) << "\n";

    std::cout << "real-1: " << rng::randomReal() << "\n";
    std::cout << "normal-1: " << rng::randomNormal() << "\n";

    RCGA ga{ 100 };
    problems::Sphere f{ 3 };

    execution_threads(1);
    rng::prng.seed(0x3da99432ab975d26LL);

    const auto solutions0 = ga.solve(f, f.bounds(), 10);
    std::cout << "single-thread: " << solutions0[0].chromosome[0] << "\n";

    execution_threads(7);
    rng::prng.seed(0x3da99432ab975d26LL);

    const auto solutions1 = ga.solve(f, f.bounds(), 10);
    std::cout << "multi-thread-1: " << solutions1[0].chromosome[0] << "\n";

    execution_threads(7);
    rng::prng.seed(0x3da99432ab975d26LL);

    const auto solutions2 = ga.solve(f, f.bounds(), 10);
    std::cout << "multi-thread-2: " << solutions2[0].chromosome[0] << "\n";
}
