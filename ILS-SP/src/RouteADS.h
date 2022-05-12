#ifndef ROUTEADS_H
#define ROUTEADS_H

#include "Data.h"
#include "infoSeq.h"
#include "Solution.h"

struct RouteADS 
{
	Data *data;
	Solution *parent;          // solution which the route belongs to 
	int cour;                  // index of the route in parent
	unsigned int length;
	double tw_violation;
	double cap_violation;
	double distance;
	double robust_tw_violation;
	double robust_cap_violation;
	double robust_demand;
	double cost;
	double min_demand;
	double max_demand;
	double min_pair_demand;
	double max_pair_demand; 
	long double random_cost;
	bool is_cap_feasible;
	bool is_time_feasible;
	bool is_feasible;
	std::vector<double> load;
	std::vector<std::pair<double, int>> sorted_deviations;
	std::vector<std::vector<Subsequence>> subseq;
	std::vector<std::vector<double>> earliest_start;
	std::vector<int> start_search;
	std::vector<int> end_search;

	RouteADS(Solution *s, int r);
	void UpdateFsr();
	void UpdateRandomCost();
	void UpdateLoadAndDistance();
	void UpdateInfoSeq();
	void UpdateRobustEarliestArrival();
	void UpdateRobustDemand();
	void UpdateMinMaxDemand();
	void Update(); // updates all ADS
};

#endif
