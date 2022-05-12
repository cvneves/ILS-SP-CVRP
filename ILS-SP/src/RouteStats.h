#ifndef ROUTESTATS_H
#define ROUTESTATS_H

#include <vector>

struct RouteStats 
{
	int nb_improvs; 
	bool is_permanent; 					
	int restart_id; 			// -1 if incumbent
	double solution_cost; // total cost of associated solution
	double solution_gap;  // gap of associated solution with respect to currently best solution 

	inline void UpdateSolutionGap(double best_solution_cost)
	{
		solution_gap = 100 * (solution_cost - best_solution_cost) / best_solution_cost;
	};

	RouteStats() {}

	RouteStats(
			int nb_improvs, 
			int restart_id, 
			bool is_permanent,
			double solution_cost, 
			double best_solution_cost
	) : nb_improvs(nb_improvs), is_permanent(is_permanent), restart_id(restart_id), solution_cost(solution_cost)
	{
		UpdateSolutionGap(best_solution_cost);
	}
};

struct RouteInfo
{
	std::vector<int> sequence;
	double cost;
	double load;
	RouteStats stats;
};

#endif
