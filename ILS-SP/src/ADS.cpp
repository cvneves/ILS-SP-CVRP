#include "ADS.h"

ADS::~ADS()
{
	for (int r = 0; r != routeADS.size(); r++)
	{
		delete routeADS[r];
	}
}

void ADS::SetSolution(Solution *s)
{
	curr_solution = s;
	for (auto r: routeADS)
	{
		delete r;	
	}

	routeADS.resize(s->route.size());
	used_vehicles = 0;

	for (int r = 0; r < s->route.size(); r++)
	{
		if (s->route[r].size() > 2)
			used_vehicles++;
		routeADS[r] = new RouteADS(s, r);
	}
}

void ADS::SetPenalizedCost()
{
	curr_solution->cost = 0;
	for (int r = 0; r < routeADS.size(); r++)
	{
		curr_solution->route_cost[r] = routeADS[r]->distance;
		curr_solution->route_load[r] = routeADS[r]->load.back();
#ifdef TW
#ifdef ROBUST_TW
		curr_solution->route_cost[r] += data->penalty_tw * routeADS[r]->robust_tw_violation;
#else
		curr_solution->route_cost[r] += data->penalty_tw * routeADS[r]->tw_violation;
#endif
#endif

#ifdef ROBUST_DEMAND
		curr_solution->route_cost[r] += data->penalty_load * routeADS[r]->robust_cap_violation;
#else
		curr_solution->route_cost[r] += data->penalty_load * routeADS[r]->cap_violation;
#endif

		curr_solution->cost += curr_solution->route_cost[r];
	}
  //   cost = 0;
  //   for (int i = 0; i < route.size(); i++)
  //   {
  //       for (int j = 1; j < route[i].size(); j++)
  //       {
  //           cost += data->time_cost[route[i][j - 1]][route[i][j]];
  //       }
  //   }
}

double ADS::GetCurrSolutionHash()
{
	double hash_val = 0;
 	std::vector<std::pair<int,int>> ord_route_idx(curr_solution->route.size());
	int j = 0;
 	for (int r = 0; r < curr_solution->route.size(); r++)
 	{
		if (curr_solution->route[r].size() <= 2)	
			continue;
// #ifndef TW	
	// 	ord_route_idx[j++] = {std::min(curr_solution->route[r][1], curr_solution->route[r].end()[-2]), r};
// #else
 		ord_route_idx[j++] = {curr_solution->route[r][1], r};
// #endif
 	}
 
 	std::sort(ord_route_idx.begin(), ord_route_idx.begin() + j);
 
 	for (int r = 0; r < j; r++)
 	{
 		hash_val += routeADS[ord_route_idx[r].second]->random_cost;
 	}
	return hash_val;
}

bool ADS::FindAndPush(int iter)
{
	double hash_val = GetCurrSolutionHash();

	auto memo_it = solution_memo.find(hash_val);
	if (memo_it == solution_memo.end()) 
	{
		solution_memo[hash_val] = iter;
		return false;
	}
	if (memo_it->second < iter)
	{ return true;	}

	return false;
}

std::pair<double, double> ADS::GetRoutePairHash(int r1, int r2)
{
	double r1_hash = routeADS[r1]->random_cost;
	double r2_hash = routeADS[r2]->random_cost;

	return {r1_hash, r2_hash};	
}

SMD* ADS::getSMD(int r1, int r2, int neighborhood)
{
	auto key = GetRoutePairHash(r1,r2);
	auto memo_it = inter_mov_memo.insert(std::pair<std::pair<double,double>, 
									std::vector<SMD>> (key, std::vector<SMD> (6, {r1,r2,0,0,0,0,0,UNVISITED,0})));
	
	return &((memo_it.first)->second[neighborhood]);
}

SMD* ADS::getSMD(int r, int neighborhood)
{
	double key = routeADS[r]->random_cost;
	auto memo_it = intra_mov_memo.insert(std::pair<double, std::vector<SMD>> (key, std::vector<SMD> (5, {r,r,0,0,0,0,0,UNVISITED,0})));

	return &((memo_it.first)->second[neighborhood]);
}
void ADS::push(int r1, int r2, int neighborhood)
{

}
void ADS::push(int r, int neighborhood)
{

}
