#ifndef STATS_H
#define STATS_H

#include "Data.h"
#include <iomanip>
#include <vector>
#include <time.h>
#include <chrono>

struct Runtime
{
	std::chrono::time_point<std::chrono::high_resolution_clock> start;
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	double total_time = 0;
	unsigned long long int count = 0;

	inline void SetStart() { start = std::chrono::high_resolution_clock::now(); }
	inline void SetEnd()
	{
		end = std::chrono::high_resolution_clock::now();

		total_time += std::chrono::duration<double>(end - start).count();
		count++;
	}
};

struct IlsStats
{
	clock_t start;
	clock_t end;
	std::vector<int> neighb_improv;
	std::vector<Runtime> inter_neighb_explore;
	std::vector<Runtime> inter_mmd_search;
	std::vector<Runtime> intra_neighb_explore;
	std::vector<Runtime> intra_mmd_search;
	std::vector<double> mmd_search_time;
	std::vector<double> inter_mmd_time;
	std::vector<double> intra_mmd_time;
	std::vector<int> inter_mmd_effective;
	std::vector<int> intra_mmd_effective;
	std::vector<int> perturb_improv;
	int last_improv_iter;

	IlsStats()
	{
		ResetNeighbImprov();
		inter_neighb_explore.assign(6, Runtime());
		intra_neighb_explore.assign(5, Runtime());
		inter_mmd_search.assign(6, Runtime());
		intra_mmd_search.assign(5, Runtime());
		inter_mmd_time.assign(6, 0);
		intra_mmd_time.assign(5, 0);
		inter_mmd_effective.assign(6, 0);
		intra_mmd_effective.assign(5, 0);
		last_improv_iter = 0;
	}
	void ResetNeighbImprov()
	{
		neighb_improv = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		perturb_improv = {0, 0};
	}
};

#endif
