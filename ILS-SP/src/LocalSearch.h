#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <bitset>
#include <fstream>
#include <functional>
#include <numeric>
#include <chrono>
#include <cmath>
#include <tuple>
#include <iomanip>

#include "Data.h"
#include "Solution.h"
#include "infoSeq.h"
#include "RouteADS.h"
#include "ADS.h"
#include "Stats.h"

#ifdef USEMMD
#define use_mmd 1
#else
#define use_mmd 0
#endif

// /** Route and Node ADS */
// struct Node
// {
//     int rank; // the node's position in the sequence
//     Node *next;
//     Node *prev;
//     Node *parent; // the route which the node belongs to
// };

// struct Route
// {
//     Data *data;
//     Node *depot;
//     unsigned int nbNodes;
//     double tw_violation;
//     double cap_violation;
//     double distance;
//     double robust_tw_violation;
//     double robust_cap_violation;
//     double robust_demand;
//     double cost;
//     bool is_cap_feasible;
//     bool is_time_feasible;
//     bool is_feasible;
//     std::vector<int> &path;
//     std::vector<double> load;
//     std::vector<std::pair<double, int>> sorted_deviations;
//     std::vector<std::vector<Subsequence>> subseq;
//     std::vector<std::vector<double>> earliest_start;

//     // Route(){};
//     Route(std::vector<int> &path) : path(path){};
//     void Update(); // updates all ADS
// };

struct ThreeBestInsertInfo
{
	//	int route[3];
	int best_position[3];
	int prev[3];
	int next[3];
	double best_cost[3] = {1e30, 1e30, 1e30};

	void CompareAndUpdateInsertion(int position, double cost, int prev, int next)
	{
		// std::cout << prev << " " << next << " " << cost << "\n";
		if (cost >= best_cost[2])
			return;
		else if (cost >= best_cost[1])
		{
			best_cost[2] = cost;
			best_position[2] = position;
			this->prev[2] = prev;
			this->next[2] = next;
		}
		else if (cost >= best_cost[0])
		{
			best_cost[2] = best_cost[1];
			best_position[2] = best_position[1];
			this->prev[2] = this->prev[1];
			this->next[2] = this->next[1];

			best_cost[1] = cost;
			best_position[1] = position;
			this->prev[1] = prev;
			this->next[1] = next;
		}
		else
		{
			best_cost[2] = best_cost[1];
			best_position[2] = best_position[1];
			this->prev[2] = this->prev[1];
			this->next[2] = this->next[1];

			best_cost[1] = best_cost[0];
			best_position[1] = best_position[0];
			this->prev[1] = this->prev[0];
			this->next[1] = this->next[0];

			best_cost[0] = cost;
			best_position[0] = position;
			this->prev[0] = prev;
			this->next[0] = next;
		}

		// std::cout << best_cost[0] << " " << best_cost[1] << " " << best_cost[2] << "\n";
	}
};

struct LocalSearch
{
	Data *data;
	ADS *ads;
	Solution *curr_solution;
	IlsStats *stats;

	std::vector<std::pair<int, int>> granSearchOrder;

	/** Local variables used during neighborhood exploration */
	std::vector<std::vector<double>> earliest_startRouteX, earliest_startRouteY, earliest_start;
	std::vector<int> pathX, pathY;

	int best_pos_x, best_pos_y;
	// int bestRouteX, bestRouteY;
	int bestMoveType;
	Subsequence subseq_x;
	Subsequence subseq_y;
	Subsequence subseq_x_rev;
	Subsequence subseq_y_rev;
	SMD *inter_smd;
	SMD *inter_smd_rev;
	SMD *intra_smd;
	int best_pair_pos_x, best_pair_pos_y;
	int move_type;
	bool improved;
	double best_delta;
	double best_rx_ry_delta;
	double curr_delta;
	double curr_delta1, curr_delta2, curr_delta3, curr_delta4;
	double temp_delta;
	double cumulated_delta;
	int cl_x, cl_z, cl_y, cl_w;
	int pos_x, pos_z, pos_y, pos_w;
	int next_z, next_w;
	int prev_x, prev_y;

	// int routeX, routeY;
	RouteADS *route_x, *route_y;	
	int rx_idx, ry_idx;
	bool rx_ry_feasible;
	bool rx_ry_improved;
	int start_search;
	int end_search;
	RouteADS *best_route_x, *best_route_y;
	double routeXRobustTWViolation, routeYRobustTWViolation;
	double routeXRobustTWViolationRev, routeYRobustTWViolationRev;
	double routeXRobustDemand, routeYRobustDemand;
	double routeXRobustCapViolation, routeYRobustCapViolation;
	int curr_neighb_idx;

	void SetLocalVariables();

	LocalSearch(Solution *s);
	LocalSearch(Solution *s, ADS *ads);

	void UpdateClientRank();

	void Run();

	void ComputeRobustDemandRouteXRelocate1();
	void ComputeRobustDemandRouteYRelocate1();
	void ConcatRouteXRelocate1Inter();
	void ConcatRouteYRelocate1Inter();
	void ConcatRouteXRelocate1IntraFront();
	void ConcatRouteXRelocate1IntraBack();
	__attribute__((always_inline)) void EvalRelocate1();
	__attribute__((always_inline)) void EvalReinsertion();
	void ExecRelocate1();
	void SearchRelocate1();
	void SearchRelocate1Complete();

	void ComputeRobustDemandRouteXRelocate2();
	void ComputeRobustDemandRouteYRelocate2();
	void ConcatRouteXRelocate2Inter();
	void ConcatRouteYRelocate2Inter();
	void ConcatRouteXRelocate2IntraFront();
	void ConcatRouteXRelocate2IntraBack();
	__attribute__((always_inline)) inline void EvalRelocate2();
	__attribute__((always_inline)) inline void EvalOrOpt();
	void ExecRelocate2();
	void SearchRelocate2();
	void SearchRelocate2Complete();

	void ComputeRobustDemandRouteXSwap11();
	void ComputeRobustDemandRouteYSwap11();
	void ConcatRouteXSwap11Inter();
	void ConcatRouteYSwap11Inter();
	void ConcatRouteXSwap11IntraFront();
	void ConcatRouteXSwap11IntraBack();
	__attribute__((always_inline)) inline void EvalSwap11();
	__attribute__((always_inline)) inline void EvalExchange();
	void ExecSwap11();
	void SearchSwap11();
	void SearchSwap11Complete();

	void ComputeRobustDemandRouteXSwap21();
	void ComputeRobustDemandRouteYSwap21();
	void ConcatRouteXSwap21Inter();
	void ConcatRouteYSwap21Inter();
	void ConcatRouteXSwap21IntraFront();
	void ConcatRouteXSwap21IntraBack();
	__attribute__((always_inline)) inline void EvalSwap21();
	void ExecSwap21();
	void SearchSwap21();
	void SearchSwap21Complete();

	void ComputeRobustDemandRouteXSwap22();
	void ComputeRobustDemandRouteYSwap22();
	void ConcatRouteXSwap22Inter();
	void ConcatRouteYSwap22Inter();
	void ConcatRouteXSwap22IntraFront();
	void ConcatRouteXSwap22IntraBack();
	__attribute__((always_inline)) inline void EvalSwap22();
	void ExecSwap22();
	void SearchSwap22();
	void SearchSwap22Complete();

	void ComputeRobustDemandRouteXCross();
	void ComputeRobustDemandRouteYCross();
	void ConcatRouteXCross();
	void ConcatRouteYCross();
	__attribute__((always_inline)) inline void EvalCross();
	void ExecCross();
	void SearchCross();
	void SearchCrossComplete();

	void ComputeRobustDemandRouteXTwoOptStar();
	void ComputeRobustDemandRouteYTwoOptStar();
	void ConcatRouteXTwoOptStar();
	void ConcatRouteYTwoOptStar();
	__attribute__((always_inline)) inline void EvalTwoOptStar();
	void ExecTwoOptStar();
	void SearchTwoOptStar();

	void ComputeRobustDemandRouteXTwoOpt();
	void ComputeRobustDemandRouteYTwoOpt();
	void ConcatRouteXTwoOpt();
	__attribute__((always_inline)) inline void EvalTwoOpt();
	void ExecTwoOpt();
	void SearchTwoOpt();

	void SearchExchange();
	void SearchReinsertion();
	void SearchOrOpt2();
	void SearchOrOpt3();
	void ExecRelocate3();

	// best_three_insertions[c][r] contains info on 2nd best insertion of client c
	// into route r
	std::vector<std::vector<ThreeBestInsertInfo>> three_best_insert; 
	int best_swapstar_ins_pos_x, best_swapstar_ins_pos_y;
	int swapstar_ins_pos_x, swapstar_ins_pos_y;
	double swapstar_ins_cost_x, swapstar_ins_cost_y;

	void ResetThreeBestInsert();
	void SearchSwapStar();
	void ComputeThreeBestInsertions();
	void ExecSwapStar();

	double ComputeCapViolation(double load)
	{
		return std::max(0.0, load - data->capacity);
	}

	void ComputePenalizedCost(Solution *s);
	void ComputeNodeEarliestStartRouteX(int currClient, int prevClient, int newRoutePos);
	void ComputeNodeEarliestStartRouteY(int currClient, int prevClient, int newRoutePos);
	void ComputePathEarliestStart(std::vector<int> &path, double &robust_tw_violation);

};

#endif
