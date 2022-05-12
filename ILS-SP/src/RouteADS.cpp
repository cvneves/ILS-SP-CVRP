#include "RouteADS.h"

RouteADS::RouteADS(Solution *s, int r)
{
	data = s->data;
	parent = s;
	cour = r;
	length = parent->route[cour].size();
	start_search.resize(data->nb_clients + 1);
	end_search.resize(data->nb_clients + 1);

	Update();
}

void RouteADS::UpdateFsr()
{
	for (int i = 1; i <= data->nb_clients; i++)
	{
		start_search[i] = 0;
		end_search[i] = length - 2;

#ifdef USEFSR
		if (is_feasible) 
		{
			for (int pos = length - 2; pos >= 1 ; pos--)
			{
				int cl = parent->route[cour][pos];
				if (!data->ins_before[i][cl])
				{
					start_search[i] = pos;
					break;
				}
			}

			for (int pos = 1; pos < length - 2; pos++)
			{
				int cl = parent->route[cour][pos];
				if (!data->ins_before[cl][i])
				{
					end_search[i] = pos - 1;
					break;
				}
			}
		}
#endif
	}
}

void RouteADS::UpdateInfoSeq()
{
	subseq = std::vector<std::vector<Subsequence>>(length, std::vector<Subsequence>(length));

	for (int i = 0; i < length; i++)
	{
		int client = parent->route[cour][i];	
		subseq[i][i].duration = data->client[client].serv_time;
		subseq[i][i].timewarp = 0;
		subseq[i][i].earliest = data->client[client].ready_time;
		subseq[i][i].latest = data->client[client].due_time;
		subseq[i][i].distance = 0;
		subseq[i][i].load = data->client[client].demand;
		subseq[i][i].first_cl = client; 
		subseq[i][i].last_cl = client;
		subseq[i][i].cap_violation = 0;

#ifdef ROBUST_CONCATENATION_TEST
		// subseq[i][i].nb_clients = path[i] != 0;
		// subseq[i][i].critLoad = std::vector<double>(data->budget_demand + 1);
		// for (int gamma = 0; gamma <= data->budget_demand; gamma++)
		// {
		// 	subseq[i][i].critLoad[gamma] = (gamma > 0) * data->client[path[i]].demand_deviation;
		// }

		// subseq[i][i].critTimewarp = std::vector<double>(data->budget_time + 1, 0.0);
		// subseq[i][i].critEarliest = std::vector<double>(data->budget_time + 1, data->client[path[i]].ready_time);
		// subseq[i][i].critLatest = std::vector<double>(data->budget_time + 1, data->client[path[i]].due_time);
		// subseq[i][i].critDuration = std::vector<double>(data->budget_time + 1, data->client[path[i]].serv_time);
#endif
	}

	for (int i = 0; i < length; i++)
		for (int j = i + 1; j < length; j++)
			subseq[i][j] = Subsequence::Concatenate(data, subseq[i][j - 1], subseq[j][j]);

	for (int i = length - 1; i >= 0; i--)
		for (int j = i - 1; j >= 0; j--)
			subseq[i][j] = Subsequence::Concatenate(data, subseq[i][j + 1], subseq[j][j]);
}

void RouteADS::UpdateRandomCost()
{
	random_cost = 0;
// #ifndef TW
// 	// CVRP instances have symmetrical graphs, so every route has an
//   // equivalent, reversed counterpart. Such route should have the same
//   // hash value, so we calculate it starting from the smallest client index.
// 	if (parent->route[cour][1] < parent->route[cour][length - 2])	
// 	{	
// 		for(int i = 0; i < length - 1; i++)
// 		{
// 			int c_i = parent->route[cour][i], c_i_next = parent->route[cour][i + 1];
// 			random_cost += data->random_cost[c_i][c_i_next];
// 		}
// 	}
// 	else
// 	{
// 		for(int i = length - 1; i > 0; i--)
// 		{
// 			int c_i = parent->route[cour][i], c_i_next = parent->route[cour][i - 1];
// 			random_cost += data->random_cost[c_i][c_i_next];
// 		}
// 	}
// #else
	for(int i = 0; i < length - 1; i++)
	{
		int c_i = parent->route[cour][i], c_i_next = parent->route[cour][i + 1];
		random_cost += (c_i + 1) * data->random_cost[c_i][c_i_next];
	}
// #endif
}

void RouteADS::UpdateLoadAndDistance()
{
	load.resize(length);

	distance = 0;
	load[0] = 0;
	for (int i = 1; i < length; i++)
	{
		load[i] = load[i-1] + data->client[parent->route[cour][i]].demand;
		distance += data->time_cost[parent->route[cour][i - 1]][parent->route[cour][i]];
	}
}

void RouteADS::UpdateRobustEarliestArrival()
{
	earliest_start = std::vector<std::vector<double>>(length, std::vector<double>(data->budget_time + 1));
	robust_tw_violation = 0;

	for (int gamma = 0; gamma <= data->budget_time; gamma++)
	{
		earliest_start[0][gamma] = 0;
	}

	for (int j = 1; j < length; j++)
	{
		int client_j = parent->route[cour][j];
		int client_j1 = parent->route[cour][j - 1];
		int currClient = parent->route[cour][j], prevClient = parent->route[cour][j - 1];
		earliest_start[j][0] = std::max(data->client[currClient].ready_time, earliest_start[j - 1][0] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient]);
	}

	for (int gamma = 1; gamma <= data->budget_time; gamma++)
	{
		for (int j = 1; j < length; j++)
		{
			int client_j = parent->route[cour][j];
			int client_j1 = parent->route[cour][j - 1];
			int currClient = parent->route[cour][j], prevClient = parent->route[cour][j - 1];
			earliest_start[j][gamma] = std::max(earliest_start[j - 1][gamma] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient], earliest_start[j - 1][gamma - 1] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient] + data->time_deviation[prevClient][currClient]);
			earliest_start[j][gamma] = std::max(earliest_start[j][gamma], data->client[currClient].ready_time);
		}
	}

	for (int j = 0; j < length; j++)
	{
		int currClient = parent->route[cour][j];
		robust_tw_violation += std::max(0.0, earliest_start[j][data->budget_time] - data->client[currClient].due_time);
	}
}

void RouteADS::UpdateRobustDemand()
{
	sorted_deviations = std::vector<std::pair<double, int>>(length);
	robust_demand = load.back();
	for (int i = 1; i < length - 1; i++)
	{
		sorted_deviations[i - 1] = {data->client[parent->route[cour][i]].demand_deviation, parent->route[cour][i]};
	}
	std::sort(sorted_deviations.begin(), sorted_deviations.end(), std::greater<std::pair<double, int>>());
	for (int i = 0; i < data->budget_demand && i < length; i++)
	{
		robust_demand += sorted_deviations[i].first;
	}
	robust_cap_violation = std::max(0.0, robust_demand - data->capacity);
}

void RouteADS::UpdateMinMaxDemand() 
{
	min_demand = std::numeric_limits<double>::infinity();
	min_pair_demand = std::numeric_limits<double>::infinity();
	max_demand = 0;
	max_pair_demand = 0;

	for (int i = 1; i < length - 1; i++)
	{
		int cl = parent->route[cour][i];
		min_demand = std::min(min_demand, data->client[cl].demand);
		max_demand = std::max(max_demand, data->client[cl].demand);
		
		if (i < length - 2)
		{
			int cl_2 = parent->route[cour][i + 1];
			min_pair_demand = std::min(min_pair_demand, data->client[cl].demand +	data->client[cl_2].demand);
			max_pair_demand = std::max(max_pair_demand, data->client[cl].demand	+ data->client[cl_2].demand);
		}
	}
}

void RouteADS::Update()
{
	/** Updating subsequence info */

	length = parent->route[cour].size();

	UpdateRandomCost();
	UpdateLoadAndDistance();
	UpdateMinMaxDemand();
	cap_violation = std::max(load.back() - data->capacity, 0.0);

#ifdef TW
	UpdateInfoSeq();

	tw_violation = subseq[0][length - 1].timewarp;

	/** robust earliest time computation */
#ifdef ROBUST_TW
	UpdateRobustEarliestArrival();
	is_time_feasible = robust_tw_violation <= EPSILON;
#else
	is_time_feasible = tw_violation <= EPSILON;
#endif
#endif

	/** Robust demands */
#ifdef ROBUST_DEMAND
	UpdateRobustDemand();
	is_cap_feasible = robust_cap_violation < EPSILON;
#else
	is_cap_feasible = cap_violation < EPSILON;
#endif

#ifdef TW
	is_feasible = (is_cap_feasible && is_time_feasible);
#else
	is_feasible = is_cap_feasible;
#endif

#ifdef TW
	UpdateFsr();
#endif
}
