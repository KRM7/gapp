/* Some test functions to compare the solutions to finding the pareto front of a set. */

#ifndef PARETO_H
#define PARETO_H

#include <vector>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <functional>
#include <execution>
#include "../src/genetic_algorithm.h"

using namespace genetic_algorithm;

std::vector<std::vector<double>> generateTestSet(size_t num_vecs, size_t dim)
{
	std::vector<std::vector<double>> ret;
	ret.reserve(num_vecs);
	for (size_t i = 0; i < num_vecs; i++)
	{
		std::vector<double> vec;
		vec.reserve(dim);
		for (size_t j = 0; j < dim; j++)
		{
			vec.push_back(rng::generateRandomDouble());
		}
		ret.push_back(vec);
	}

	return ret;
}

bool paretoCompare(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
	bool lhs_is_dominated = false;
	for (size_t i = 0; i < lhs.size(); i++)
	{
		if (lhs[i] > rhs[i]) return false;
		if (lhs[i] < rhs[i]) lhs_is_dominated = true;
	}

	return lhs_is_dominated;
}

std::vector<std::vector<double>> findParetoFrontKung(std::vector<std::vector<double>>& pop)
{
	/* Kung's algorithm. */
	/* See: Kung et al. "On finding the maxima of a set of vectors." Journal of the ACM (JACM) 22.4 (1975): 469-476. */	
	using namespace std;
	using iter = vector<size_t>::iterator;

	/* Sort pop into descending order based on first element. */
	sort(pop.begin(), pop.end(), [](const vector<double>& lhs, const vector<double>& rhs) {return lhs[0] > rhs[0]; });

	/* Finds the indices of pareto optimal vectors in pop. */
	function<vector<size_t>(iter, iter)> pfront = [&pop, &pfront](iter first, iter last) -> vector<size_t>
	{
		if (distance(first, last) == 1) return { *first };

		vector<size_t> R = pfront(first, first + distance(first, last) / 2);	/* Top half. */
		vector<size_t> S = pfront(first + distance(first, last) / 2, last);		/* Bottom half. */

		/* T = Find all s elements of S which are not dominated by any r element of R. */
		vector<size_t> T;
		for (size_t s = 0; s < S.size(); s++)
		{
			bool is_dominated = false;
			for (size_t r = 0; r < R.size(); r++)
			{
				/* Pareto comparison between s and r. */
				/* The first dimension of the vectors doesn't need to be compared since it's already sorted. */
				for (size_t k = 1; k < pop[s].size(); k++)
				{
					if (pop[S[s]][k] > pop[R[r]][k])
					{
						is_dominated = false;
						break;	/* If s has a strictly better element, it can't be dominated by r. */
					}
					if (pop[S[s]][k] < pop[R[r]][k]) is_dominated = true;
				}
				if (is_dominated) break;	/* Already dominated, no reason to check the rest. */
			}
			if (!is_dominated) T.push_back(S[s]);
		}
		R.insert(R.end(), T.begin(), T.end());

		return R;
	};

	vector<size_t> indices(pop.size());
	iota(indices.begin(), indices.end(), 0U);

	indices = pfront(indices.begin(), indices.end());

	vector<vector<double>> pareto_sols;
	pareto_sols.reserve(indices.size());
	for (size_t i = 0; i < indices.size(); i++)
	{
		pareto_sols.push_back(pop[indices[i]]);
	}

	return pareto_sols;
}

std::vector<std::vector<double>> findParetoFrontNaive(const std::vector<std::vector<double>>& pop)
{
	std::vector<std::vector<double>> optimal_sols;
	optimal_sols.reserve(pop.size());

	for (auto lhs = pop.begin(); lhs != pop.end(); lhs++)
	{
		bool dominated = false;
		for (auto rhs = pop.begin(); rhs != pop.end(); rhs++)
		{
			dominated = paretoCompare(*lhs, *rhs);
			if (dominated) break;
		}
		if (!dominated) optimal_sols.push_back(*lhs);
	}

	return optimal_sols;
}

/* Return a vector of indices for each front. */
std::vector<std::vector<size_t>> findParetoFronts(const std::vector<std::vector<double>>& pop)
{
	using namespace std;

	/* Calc the number of vectors which dominate each vector (dom_count), and the set of dominated vectors for each vector (dom_list). */
	vector<size_t> dom_count(pop.size(), 0U);
	vector<vector<size_t>> dom_list(pop.size());

	for (size_t i = 0; i < pop.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (paretoCompare(pop[j], pop[i]))	/* i dominates j */
			{
				dom_count[j]++;
				dom_list[i].push_back(j);
			}
			else if (paretoCompare(pop[i], pop[j])) /* j dominates i */
			{
				dom_count[i]++;
				dom_list[j].push_back(i);
			}
		}
	}

	/* Find all vectors which are not dominated by any vectors (first pareto front). */
	vector<size_t> front;
	for (size_t i = 0; i < pop.size(); i++)
	{
		if (dom_count[i] == 0)
		{
			front.push_back(i);
			//assign rank 0 to pop[i]
		}
	}

	/* Find all the other fronts. */
	vector<vector<size_t>> pareto_fronts;
	size_t front_idx = 1U;
	while (!front.empty())
	{
		/* "Remove" the current front and find the next one. */
		vector<size_t> next_front;
		for (auto i = front.begin(); i != front.end(); i++)
		{
			for (auto j = dom_list[*i].begin(); j != dom_list[*i].end(); j++)	/* Decrease dom_count of every vector dominated by i. */
			{
				dom_count[*j]--;
				if (dom_count[*j] == 0)	/* j belongs to the next front if its domination count became 0. */
				{
					next_front.push_back(*j);
					//assign rank front_idx to pop[*j]
				}
			}
		}
		pareto_fronts.push_back(front);
		front = next_front;
		front_idx++;
	}

	return pareto_fronts;
}

void testParetoFront(size_t dim)
{
	using namespace std;
	using namespace chrono;

	for (size_t n = 100; n <= 12800; n *= 2)
	{
		vector<vector<double>> test_set = generateTestSet(n, dim);

		cout << fixed << setprecision(4);

		/* Naive algorithm. */
		auto t_begin = high_resolution_clock::now();
		vector<std::vector<double>> res1 = findParetoFrontNaive(test_set);
		auto t_end = high_resolution_clock::now();

		auto duration = duration_cast<milliseconds>(t_end - t_begin).count();
		double ts = duration / 1E+3;
		cout << "Naive algorithm for " << n << " elements in " << dim << " dimensions:\t" << ts << " s\n";

		/* Kung's algorithm. */
		t_begin = high_resolution_clock::now();
		vector<std::vector<double>> res2 = findParetoFrontKung(test_set);
		t_end = high_resolution_clock::now();

		duration = duration_cast<milliseconds>(t_end - t_begin).count();
		ts = duration / 1E+3;
		cout << "Kung's algorithm for " << n << " elements in " << dim << " dimensions:\t" << ts << " s\n";

		cout << "Same results: " << is_permutation(res1.begin(), res1.end(), res2.begin()) << "\n\n";

		/*
		* The naive algorithm is faster for d = 2 and every reasonable n.
		* Kung's is faster for:
		*	d = 3 and n > 1600.
		*	d = 4 and n > 1600.
		*	d = 6 and n > 800.
		*	d = 10 and n > 400.
		* Kung's is faster for every n in large dimensions (d > 100).
		* 
		* The main advantage of Kung's algorithm is for large numbers of solutions (n > 1000) or large dimensions (d > 20),
		* otherwise it has around the same performance as the naive algorithm.
		*/
	}
}

#endif // !PARETO_H