#ifndef ILSRVND_H
#define ILSRVND_H

#include <iostream>
#include <fstream>
#include "Data.h"
#include "Solution.h"
#include "LocalSearch.h"
#include "Pool.h"
#include "SP.h"
#include "Construction.h"
#include "Stats.h"
#include "ADS.h"
#include "Perturbation.h"

#include <string>
#include <map>
#include <deque>
#include <tuple>
#include <chrono>

using namespace std;
struct SPSolver;

struct AdaptativeInfo
{
		Data *data;

		int reference_sol_idx = 0;
		int last_sol_idx = 0;
		int best_sol_idx = 0;
		size_t max_nb_solutions;
		double threshold = std::numeric_limits<double>::infinity();
		double quality_factor = 1.0;
		double decrease_factor;
		double average_cost = 0;

		deque<Solution> previous_solutions;

		void Update(Solution *solution);
		void UpdatePrevSolutions(Solution *solution);
		void UpdateBestSolutionIdx();
		void UpdateQualityFactor();
		void UpdateThreshold(Solution *solution);

		AdaptativeInfo(size_t max_solutions, double decrease_factor);
		Solution GetSolution();
};

struct IlsRvnd
{
	Data *data;
	Pool *my_pool = NULL;
	SPSolver *my_sp;

	Solution best;
	Solution *init_solution;

	bool skip_construction = false;

	int iteration = 0;
	int ils_improv;
	double tolerance;
	double init_sol_cost;
	static double best_of_all_cost;
	static int iter_id;

	double ils_gap;
	double sp_gap;
	
	std::deque<std::tuple<double, double>> cost_history;
	size_t cost_history_size = 1250;
	size_t cost_history_step = 500;
	size_t history_counter;

	GlobalSolution *global_solution = NULL;
	int *clust_idx = NULL;

	IlsRvnd(Data *data);
	IlsRvnd(Data *data, Pool *my_pool, Solution *s0);
	IlsRvnd(Data *data, Pool *my_pool);

	void Run(int maxIter, int max_iter_ils, double tolerance);
	void LogCurrentIterationStats(IlsStats *stats);
	void LogFinalStats(IlsStats *stats);

	void Improve();

	int EstimateNumberOfVehicles();
};

#endif
