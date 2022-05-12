#include "LocalSearch.h"

void LocalSearch::ComputePenalizedCost(Solution *solution)
{
	//     solution->cost = 0;
	//     for (int r = 0; r < ads->routeADS.size(); r++)
	//     {
	//         solution->route_cost[r] = routeData[r].distance;
	// #ifdef TW
	// #ifdef ROBUST_TW
	//         solution->route_cost[r] += data->penalty_tw * routeData[r].robust_tw_violation;
	// #else
	//         solution->route_cost[r] += data->penalty_tw * routeData[r].tw_violation;
	// #endif
	// #endif
	// 
	// #ifdef ROBUST_DEMAND
	//         solution->route_cost[r] += data->penalty_load * routeData[r].robust_cap_violation;
	// #else
	//         solution->route_cost[r] += data->penalty_load * routeData[r].cap_violation;
	// #endif
	// 
	//         solution->cost += solution->route_cost[r];
	//     }
}

void LocalSearch::UpdateClientRank()
{
	// 	for (int r = 0; r < ads->routeADS.size(); r++)
	// 		for (int i = 1; i < routeData[r].path.size() - 1; i++)
	// 			clientRank[routeData[r].path[i]] = {r, i};
}

LocalSearch::LocalSearch(Solution *solution) : data(solution->data), curr_solution(solution)
{
	//     routeData.clear();
	//     clientRank = std::vector<std::pair<int, int>>(data->nb_clients + 1);
	// 
	//     for (int r = 0; r < curr_solution->route.size(); r++)
	//     {
	//         routeData.push_back(Route(curr_solution->route[r]));
	//         routeData[r].data = curr_solution->data;
	//         // routeData[r].path = solution->route[r];
	//         routeData[r].Update();
	//     }
	//     UpdateClientRank();
	// 
	//     for (int i = 0; i < data->correlated_vertices.size(); i++)
	//         for (const auto &j : data->correlated_vertices[i])
	//             granSearchOrder.push_back({i, j});
	// 
	//     earliest_startRouteX = std::vector<std::vector<double>>(data->nb_clients + 2, std::vector<double>(data->budget_time + 1));
	//     earliest_startRouteY = std::vector<std::vector<double>>(data->nb_clients + 2, std::vector<double>(data->budget_time + 1));
	//     earliest_start = std::vector<std::vector<double>>(data->nb_clients + 2, std::vector<double>(data->budget_time + 1));
	// 
	//     // ComputePenalizedCost(curr_solution);
	//     // for(int i = 0; i < ads->routeADS.size(); i++)
	//     // {
	//     //     std::cout << routeData[i].subseq[0][routeData[i].path.size() - 1].distance << "\n";
	//     // }
}

LocalSearch::LocalSearch(Solution *solution, ADS *ads) : ads(ads), data(solution->data), curr_solution(solution), stats(NULL)
{
#ifdef ROBUST_TW
	earliest_startRouteX = std::vector<std::vector<double>>(data->nb_clients + 2, std::vector<double>(data->budget_time + 1));
	earliest_startRouteY = std::vector<std::vector<double>>(data->nb_clients + 2, std::vector<double>(data->budget_time + 1));
	earliest_start       = std::vector<std::vector<double>>(data->nb_clients + 2, std::vector<double>(data->budget_time + 1));
#endif

	// ComputePenalizedCost(curr_solution);
	// for(int i = 0; i < ads->routeADS.size(); i++)
	// {
	//     std::cout << routeData[i].subseq[0][routeData[i].path.size() - 1].distance << "\n";
	// }
}

void LocalSearch::Run()
{
	// std::cout << "START\n";
	// std::vector<void (LocalSearch::*)()> NL = {&LocalSearch::SearchRelocate1Complete, &LocalSearch::SearchRelocate2Complete, &LocalSearch::SearchSwap11Complete, &LocalSearch::SearchSwap21Complete, &LocalSearch::SearchSwap22Complete, &LocalSearch::SearchCrossComplete, &LocalSearch::SearchSwapStar};

	//
	// ads->inter_mov_memo.clear();
	// ads->intra_mov_memo.clear();

	//std::vector<void (LocalSearch::*)()> NL = {&LocalSearch::SearchCrossComplete};
	std::vector<void (LocalSearch::*)()> NL = {&LocalSearch::SearchRelocate1Complete,
	&LocalSearch::SearchRelocate2Complete,
	&LocalSearch::SearchSwap11Complete,
	&LocalSearch::SearchSwap22Complete,
	&LocalSearch::SearchSwap21Complete,
	&LocalSearch::SearchCrossComplete};

	//&LocalSearch::SearchSwapStar};

	if (data->use_swap_star)
		NL.push_back(&LocalSearch::SearchSwapStar);

	//std::vector<void (LocalSearch::*)()> intraNL = {&LocalSearch::SearchTwoOpt};
	std::vector<void (LocalSearch::*)()> intraNL = {&LocalSearch::SearchExchange, 
	&LocalSearch::SearchTwoOpt,
	&LocalSearch::SearchReinsertion,
	&LocalSearch::SearchOrOpt2,
	&LocalSearch::SearchOrOpt3};

	std::vector<void (LocalSearch::*)()> curr_nl = NL;

	int iter = 0;
	cumulated_delta = 0;

	int conta = 0;
	while (!curr_nl.empty())
	{
		best_route_x = best_route_y = NULL;
		conta++;
		int neighb_idx = rand() % curr_nl.size();
		// bestMoveType = 10;

		(this->*curr_nl[neighb_idx])();

		if (improved)
		{
			curr_solution->SetIntraNeighbStatusAllTrue(best_route_x->cour);
			curr_solution->SetIntraNeighbStatusAllTrue(best_route_y->cour);
			cumulated_delta += best_delta;
			curr_nl = NL;

			std::vector<void (LocalSearch::*)()> curr_intra_nl = intraNL;
			
			while (!curr_intra_nl.empty())
			{
				int intra_neighb_idx = rand() % curr_intra_nl.size();
				(this->*curr_intra_nl[intra_neighb_idx])();

				curr_intra_nl.erase(curr_intra_nl.begin() + intra_neighb_idx);
			}
		}
		else
		{
			curr_nl.erase(curr_nl.begin() + neighb_idx);
		}
	}
}

inline void LocalSearch::SetLocalVariables()
{
	//     pos_z = pos_x = clientRank[cl_x].second;
	//     pos_w = pos_y = clientRank[cl_y].second;
	//     pos_z++, pos_w++;
	//     routeX = clientRank[cl_x].first;
	//     routeY = clientRank[cl_y].first;
	//     cl_z = curr_solution->route[rx_idx][pos_z];
	//     cl_w = curr_solution->route[ry_idx][pos_w];
	//     prev_x = curr_solution->route[rx_idx][pos_x - 1];
	//     prev_y = curr_solution->route[ry_idx][pos_y - (pos_y > 0)];
}

inline void LocalSearch::ComputeNodeEarliestStartRouteX(int currClient, int prevClient, int newRoutePos)
{
	if (newRoutePos == 0)
	{
		for (int gamma = 0; gamma <= data->budget_time; gamma++)
		{
			earliest_startRouteX[newRoutePos][gamma] = 0;
		}
		routeXRobustTWViolation = 0;
	}
	else
	{
		earliest_startRouteX[newRoutePos][0] = std::max(data->client[currClient].ready_time, earliest_startRouteX[newRoutePos - 1][0] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient]);
		for (int gamma = 1; gamma <= data->budget_time; gamma++)
		{
			earliest_startRouteX[newRoutePos][gamma] = std::max(data->client[currClient].ready_time, earliest_startRouteX[newRoutePos - 1][gamma] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient]);
			earliest_startRouteX[newRoutePos][gamma] = std::max(earliest_startRouteX[newRoutePos][gamma], earliest_startRouteX[newRoutePos - 1][gamma - 1] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient] + data->time_deviation[prevClient][currClient]);
		}
		routeXRobustTWViolation += std::max(0.0, earliest_startRouteX[newRoutePos][data->budget_time] - data->client[currClient].due_time);
	}
}

inline void LocalSearch::ComputeNodeEarliestStartRouteY(int currClient, int prevClient, int newRoutePos)
{
	if (newRoutePos == 0)
	{
		for (int gamma = 0; gamma <= data->budget_time; gamma++)
		{
			earliest_startRouteY[newRoutePos][gamma] = 0;
		}
		routeYRobustTWViolation = 0;
	}
	else
	{
		earliest_startRouteX[newRoutePos][0] = std::max(data->client[currClient].ready_time, earliest_startRouteY[newRoutePos - 1][0] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient]);
		for (int gamma = 1; gamma <= data->budget_time; gamma++)
		{
			earliest_startRouteY[newRoutePos][gamma] = std::max(data->client[currClient].ready_time, earliest_startRouteY[newRoutePos - 1][gamma] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient]);
			earliest_startRouteY[newRoutePos][gamma] = std::max(earliest_startRouteY[newRoutePos][gamma], earliest_startRouteY[newRoutePos - 1][gamma - 1] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient] + data->time_deviation[prevClient][currClient]);
		}
		routeYRobustTWViolation += std::max(0.0, earliest_startRouteY[newRoutePos][data->budget_time] - data->client[currClient].due_time);
	}
}

inline void LocalSearch::ComputePathEarliestStart(std::vector<int> &path, double &robust_tw_violation)
{
	int rLen = path.size();
	for (int gamma = 0; gamma <= data->budget_time; gamma++)
	{
		earliest_start[0][gamma] = 0;
	}

	for (int j = 1; j < rLen; j++)
	{
		int currClient = path[j], prevClient = path[j - 1];
		earliest_start[j][0] = std::max(data->client[currClient].ready_time, earliest_start[j - 1][0] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient]);
	}

	for (int gamma = 1; gamma <= data->budget_time; gamma++)
	{
		for (int j = 1; j < rLen; j++)
		{
			int currClient = path[j], prevClient = path[j - 1];
			earliest_start[j][gamma] = std::max(earliest_start[j - 1][gamma] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient], earliest_start[j - 1][gamma - 1] + data->client[prevClient].serv_time + data->time_cost[prevClient][currClient] + data->time_deviation[prevClient][currClient]);
			earliest_start[j][gamma] = std::max(earliest_start[j][gamma], data->client[currClient].ready_time);
		}
	}

	for (int j = 0; j < rLen; j++)
	{
		int currClient = path[j];
		robust_tw_violation += std::max(0.0, earliest_start[j][data->budget_time] - data->client[currClient].due_time);
	}
}

/** ********************* RELOCATE 1 ***************************** */

inline void LocalSearch::ComputeRobustDemandRouteXRelocate1()
{
	routeXRobustDemand = route_x->load.back() - data->client[cl_x].demand;
	std::vector<double> sorted_deviations(route_x->length - 3);

	for (int i = 1, j = 0; i <= route_x->length - 2; i++)
	{
		if (i == pos_x)
			continue;
		sorted_deviations[j] = -data->client[curr_solution->route[rx_idx][i]].demand_deviation;
		j++;
	}

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}

	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}

	routeXRobustCapViolation = std::max(0.0, routeXRobustDemand - data->capacity);
}

inline void LocalSearch::ComputeRobustDemandRouteYRelocate1()
{
	double routeYRobustDemand = route_y->load.back() + data->client[cl_x].demand;
	std::vector<double> sorted_deviations(route_y->length - 1);

	int i;
	for (i = 1; i <= route_y->length - 2; i++)
		sorted_deviations[i - 1] = -data->client[curr_solution->route[ry_idx][i]].demand_deviation;
	sorted_deviations.back() = -data->client[cl_x].demand_deviation;

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}

	routeYRobustCapViolation = std::max(0.0, routeYRobustDemand - data->capacity);
}

inline void LocalSearch::ConcatRouteXRelocate1Inter()
{
#ifdef ROBUST_TW
	// int l = 0;
	// ComputeNodeEarliestStartRouteX(curr_solution->route[rx_idx][0], 0, l++);

	// for (int j = 1; j <= pos_x - 1; j++)
	//     ComputeNodeEarliestStartRouteX(curr_solution->route[rx_idx][j], curr_solution->route[rx_idx][j - 1], l++);

	// ComputeNodeEarliestStartRouteX(curr_solution->route[rx_idx][pos_x + 1], curr_solution->route[rx_idx][pos_x - 1], l++);

	// for (int j = pos_x + 2; j <= route_x->length - 1; j++)
	//     ComputeNodeEarliestStartRouteX(curr_solution->route[rx_idx][j], curr_solution->route[rx_idx][j - 1], l++);

	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_x + 1; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_x + 1][route_x->length - 1]);
}

inline void LocalSearch::ConcatRouteYRelocate1Inter()
{
#ifdef ROBUST_TW
	// int l = 0;
	// ComputeNodeEarliestStartRouteY(curr_solution->route[ry_idx][0], 0, l++);
	// for (int j = 1; j <= pos_y; j++)
	//     ComputeNodeEarliestStartRouteY(curr_solution->route[ry_idx][j], curr_solution->route[ry_idx][j - 1], l++);

	// ComputeNodeEarliestStartRouteY(curr_solution->route[rx_idx][pos_x], curr_solution->route[ry_idx][pos_y], l++);
	// ComputeNodeEarliestStartRouteY(curr_solution->route[ry_idx][pos_y + 1], curr_solution->route[rx_idx][pos_x], l++);

	// for (int j = pos_y + 2; j <= route_y->length - 1; j++)
	//     ComputeNodeEarliestStartRouteY(curr_solution->route[ry_idx][j], curr_solution->route[ry_idx][j - 1], l++);
	routeYRobustTWViolation = 0;
	pathY.clear();
	for (int i = 0; i <= pos_y; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	pathY.push_back(curr_solution->route[rx_idx][pos_x]);
	for (int i = pos_y + 1; i < route_y->length; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	ComputePathEarliestStart(pathY, routeYRobustTWViolation);
#endif

	subseq_y = Subsequence::Concatenate(data, route_y->subseq[0][pos_y], route_x->subseq[pos_x][pos_x]);
	subseq_y = Subsequence::Concatenate(data, subseq_y, route_y->subseq[pos_y + 1][route_y->length - 1]);
}

inline void LocalSearch::ConcatRouteXRelocate1IntraFront()
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_x + 1; i <= pos_y; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_x);
	for (int i = pos_y + 1; i <= route_x->length - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}

	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_x + 1][pos_y]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x][pos_x]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][route_x->length - 1]);
}

inline void LocalSearch::ConcatRouteXRelocate1IntraBack()
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_y; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_x);
	for (int i = pos_y + 1; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_x + 1; i <= route_x->length - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}

	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_y], route_x->subseq[pos_x][pos_x]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][pos_x - 1]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 1][route_x->length - 1]);
}

inline void LocalSearch::EvalRelocate1()
{
	curr_delta = temp_delta - data->time_cost[cl_y][cl_w] + data->time_cost[cl_y][cl_x] + data->time_cost[cl_x][cl_w];

	// if (route_x->is_feasible && route_y->is_feasible && curr_delta >= EPSILON)
	// 	return;

#ifdef TW
	ConcatRouteXRelocate1Inter();
	ConcatRouteYRelocate1Inter();

#ifdef ROBUST_TW
	curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolation);
#else
	curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp) +
		data->penalty_tw * (-route_y->tw_violation + subseq_y.timewarp);
#endif
#endif

#ifdef ROBUST_DEMAND
	ComputeRobustDemandRouteXRelocate1();
	ComputeRobustDemandRouteYRelocate1();

	curr_delta += data->penalty_load * (routeXRobustCapViolation + routeYRobustCapViolation -
			route_x->robust_cap_violation - route_y->robust_cap_violation);
#else
	//             curr_delta += data->penalty_load * (-route_x->cap_violation + subseq_x.cap_violation) +
	//                          data->penalty_load * (-route_y->cap_violation + subseq_y.cap_violation);

	curr_delta += data->penalty_load * (-route_x->cap_violation + ComputeCapViolation(route_x->load.back() - data->client[cl_x].demand)) +
		data->penalty_load * (-route_y->cap_violation + ComputeCapViolation(route_y->load.back() + data->client[cl_x].demand));
#endif

#ifdef USEMMD
	if (curr_delta + EPSILON < inter_smd->cost)
	{
		inter_smd->status = IMPROVED;
		rx_ry_improved = true;
		inter_smd->cost = curr_delta;
		inter_smd->r1_pos1 = pos_x;
		inter_smd->r2_pos1 = pos_y;
	}
#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
	}
#endif
}

inline void LocalSearch::EvalReinsertion()
{
	curr_delta = -data->time_cost[prev_x][cl_x] - data->time_cost[cl_x][cl_z] + data->time_cost[prev_x][cl_z] -
		data->time_cost[cl_y][cl_w] + data->time_cost[cl_y][cl_x] + data->time_cost[cl_x][cl_w];
#ifdef TW
	if (route_x->is_time_feasible && curr_delta >= EPSILON)
		return;

	if (pos_x < pos_y)
		ConcatRouteXRelocate1IntraFront();
	else
		ConcatRouteXRelocate1IntraBack();

#ifdef ROBUST_TW
	curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation);
#else
	curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
#endif
#endif

#ifdef USEMMD
	if (curr_delta + EPSILON < intra_smd->cost)
	{
		intra_smd->status = IMPROVED;
		intra_smd->cost = curr_delta;
		intra_smd->r1_pos1 = pos_x;
		intra_smd->r2_pos1 = pos_y;
	}
#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
	}
#endif
}

inline void LocalSearch::ExecRelocate1()
{
	if (best_route_x != best_route_y)
	{
		for (int i = best_pos_x; i < best_route_x->length - 1; i++)
			std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);
		int cust = curr_solution->route[best_route_x->cour].back();
		curr_solution->route[best_route_x->cour].pop_back();
		best_route_x->Update();

		curr_solution->route[best_route_y->cour].push_back(cust);
		for (int i = best_route_y->length - 1; i > best_pos_y; i--)
			std::swap(curr_solution->route[best_route_y->cour][i], curr_solution->route[best_route_y->cour][i + 1]);
		best_route_y->Update();
	}
	else
	{
		if (best_pos_x < best_pos_y)
			for (int i = best_pos_x; i < best_pos_y; i++)
				std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);
		else
			for (int i = best_pos_x - 1; i > best_pos_y; i--)
				std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);
		best_route_x->Update();
	}

	// UpdateClientRank();
}

void LocalSearch::SearchRelocate1()
{
	//     improved = false;
	//     best_delta = 0;
	//     for (cl_x = 1; cl_x <= data->nb_clients; cl_x++)
	//     {
	//         for (int cy = 0; cy < data->correlated_vertices[cl_x].size(); cy++)
	//         {
	//             cl_y = data->correlated_vertices[cl_x][cy];
	//             SetLocalVariables();
	//             EvalRelocate1();
	//             // if (pos_y == 1)
	//             // {
	//             //     pos_y--;
	//             //     pos_w--;
	// 
	//             //     cl_y = 0;
	//             //     cl_w = curr_solution->route[ry_idx][pos_w];
	//             //     prev_y = curr_solution->route[ry_idx][pos_y - (pos_y > 0)];
	//             //     EvalRelocate1();
	//             // }
	//             // std::cout << *curr_solution << "\n";
	//         }
	//     }
	// 
	//     if (improved)
	//     {
	//         // std::cout << *curr_solution << "\n";
	//         // std::cout << best_delta << "\n";
	//         // std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";
	// 
	//         ExecRelocate1();
	//     }
}

void LocalSearch::SearchRelocate1Complete()
{
	stats->inter_neighb_explore[0].SetStart();

	curr_neighb_idx = 0;
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		route_x = ads->routeADS[rx_idx];
		if (route_x->length <= 2)
			continue;

		for (ry_idx = 0; ry_idx < ads->routeADS.size(); ry_idx++)
		{
			route_y = ads->routeADS[ry_idx];
			rx_ry_feasible = route_x->is_feasible && route_y->is_feasible;

			if (rx_idx == ry_idx
				|| (rx_ry_feasible && (route_y->load.back() + route_x->min_demand > data->capacity))
				|| (!curr_solution->GetInterNeighbStatus(rx_idx, ry_idx, 0)))
				continue;
			
			rx_ry_improved  = false;

			stats->inter_mmd_search[0].SetStart();
#ifdef USEMMD
			inter_smd = ads->getSMD(rx_idx, ry_idx, 0);	
#endif
			stats->inter_mmd_search[0].SetEnd();

			start_search = 0;
			end_search  = route_y->length - 2;

			if (!use_mmd || inter_smd->status == UNVISITED)
			{
#ifdef USEMMD
				inter_smd->status = UNIMPROVED;
#endif


				for (pos_x = 1; pos_x < route_x->length - 1; pos_x++)
				{
					cl_x = curr_solution->route[rx_idx][pos_x];
					if (route_y->length <= 2 || rx_ry_feasible && (route_y->load.back() + data->client[cl_x].demand > data->capacity))
						continue;

					pos_z = pos_x + 1;
					cl_z = curr_solution->route[rx_idx][pos_z];
					prev_x = curr_solution->route[rx_idx][pos_x - 1];

					temp_delta = - data->time_cost[prev_x][cl_x] - data->time_cost[cl_x][cl_z] + data->time_cost[prev_x][cl_z];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_x];
						end_search = route_y->end_search[cl_x];
					}
#endif

					for (pos_y = start_search; pos_y <= end_search; pos_y++)
					{
						pos_w = pos_y + 1;
						cl_y = curr_solution->route[ry_idx][pos_y];
						cl_w = curr_solution->route[ry_idx][pos_w];
						// prev_y = curr_solution->route[ry_idx][pos_y - (pos_y > 0)];

						EvalRelocate1();
						rx_ry_improved |= (curr_delta + EPSILON < 0) ;
					}
				}

			}		
			else
			{
				stats->inter_mmd_effective[0]++;
			}

#ifdef USEMMD		
			if (inter_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = inter_smd->cost;
				best_pos_x = inter_smd->r1_pos1;
				best_pos_y = inter_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
			}
			rx_ry_improved |= (inter_smd->cost + EPSILON < 0);
#endif

			if (!rx_ry_improved)
			 curr_solution->SetInterNeighbStatus(rx_idx, ry_idx, 0, false);
		}
	}

	if (improved)
	{
		stats->neighb_improv[0]++;
		 // std::cout << *curr_solution << "\n";
		 // std::cout << "Best delta: " << best_delta << "\n";
		 // std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x->cour << " " << best_route_y->cour << "\n";

		curr_solution->SetInterNeighbStatusAllTrue(best_route_x->cour, best_route_y->cour);
		ExecRelocate1();
	}
	stats->inter_neighb_explore[0].SetEnd();
}

/** ********************* SWAP (1,1) ***************************** */

inline void LocalSearch::ComputeRobustDemandRouteXSwap11()
{
	routeXRobustDemand = route_x->load.back() - data->client[cl_x].demand + data->client[cl_y].demand;
	std::vector<double> sorted_deviations(route_x->length - 2);

	int i;
	for (i = 1; i <= route_x->length - 2; i++)
	{
		if (i == pos_x)
		{
			sorted_deviations[i - 1] = -data->client[cl_y].demand;
			continue;
		}
		sorted_deviations[i - 1] = -data->client[curr_solution->route[rx_idx][i]].demand_deviation;
	}

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}
	routeXRobustCapViolation = std::max(0.0, routeXRobustDemand - data->capacity);
}

inline void LocalSearch::ComputeRobustDemandRouteYSwap11()
{
	double routeYRobustDemand = route_y->load.back() - data->client[cl_y].demand + data->client[cl_x].demand;
	std::vector<double> sorted_deviations(route_y->length - 2);

	int i;
	for (i = 1; i <= route_y->length - 2; i++)
	{
		if (i == pos_y)
		{
			sorted_deviations[i - 1] = -data->client[cl_x].demand;
			continue;
		}
		sorted_deviations[i - 1] = -data->client[curr_solution->route[ry_idx][i]].demand_deviation;
	}

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}

	routeYRobustCapViolation = std::max(0.0, routeYRobustDemand - data->capacity);
}

inline void LocalSearch::ConcatRouteXSwap11Inter()
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	for (int i = pos_x + 1; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_y->subseq[pos_y][pos_y]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 1][route_x->length - 1]);
}

inline void LocalSearch::ConcatRouteYSwap11Inter()
{
#ifdef ROBUST_TW
	routeYRobustTWViolation = 0;
	pathY.clear();
	for (int i = 0; i <= pos_y - 1; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	pathY.push_back(cl_x);
	for (int i = pos_y + 1; i < route_y->length; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	ComputePathEarliestStart(pathY, routeYRobustTWViolation);
#endif

	subseq_y = Subsequence::Concatenate(data, route_y->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x]);
	subseq_y = Subsequence::Concatenate(data, subseq_y, route_y->subseq[pos_y + 1][route_y->length - 1]);
}

inline void LocalSearch::ConcatRouteXSwap11IntraFront() // when pos_x < pos_y
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	for (int i = pos_x + 1; i <= pos_y - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_x);
	for (int i = pos_y + 1; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif
	if (pos_x + 1 == pos_y)
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_y][pos_y]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x][pos_x]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][route_x->length - 1]);
	}
	else
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_y][pos_y]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 1][pos_y - 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x][pos_x]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][route_x->length - 1]);
	}
}

inline void LocalSearch::ConcatRouteXSwap11IntraBack() // when pos_x > pos_y
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_y - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_x);
	for (int i = pos_y + 1; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	for (int i = pos_x + 1; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	if (pos_y + 1 == pos_x)
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y][pos_y]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 1][route_x->length - 1]);
	}
	else
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][pos_x - 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y][pos_y]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 1][route_x->length - 1]);
	}
}

inline void LocalSearch::EvalSwap11()
{
	curr_delta = -data->time_cost[prev_x][cl_x] - data->time_cost[cl_x][cl_z] + data->time_cost[prev_x][cl_y] + data->time_cost[cl_y][cl_z] +
		-data->time_cost[prev_y][cl_y] - data->time_cost[cl_y][cl_w] + data->time_cost[prev_y][cl_x] + data->time_cost[cl_x][cl_w];

	if (route_x->is_feasible && route_y->is_feasible && curr_delta >= EPSILON)
		return;

#ifdef TW
	ConcatRouteXSwap11Inter();
	ConcatRouteYSwap11Inter();

#ifdef ROBUST_TW
	curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolation);
#else
	curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp) +
		data->penalty_tw * (-route_y->tw_violation + subseq_y.timewarp);
#endif
#endif

#ifdef ROBUST_DEMAND
	ComputeRobustDemandRouteXSwap11();
	ComputeRobustDemandRouteYSwap11();
	curr_delta += data->penalty_load * (routeXRobustCapViolation + routeYRobustCapViolation -
			route_x->robust_cap_violation - route_y->robust_cap_violation);
#else
	double tempLoadDelta = -data->client[cl_x].demand + data->client[cl_y].demand;
	curr_delta += data->penalty_load * (-route_x->cap_violation + ComputeCapViolation(route_x->load.back() + tempLoadDelta)) +
		data->penalty_load * (-route_y->cap_violation + ComputeCapViolation(route_y->load.back() - tempLoadDelta));
#endif

#ifdef USEMMD
	if (curr_delta + EPSILON < inter_smd->cost)
	{
		inter_smd->status = IMPROVED;
		inter_smd->cost = curr_delta;
		inter_smd->r1_pos1 = pos_x;
		inter_smd->r2_pos1 = pos_y;
	}
#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
	}
#endif
}

inline void LocalSearch::EvalExchange()
{
	// std::cout << pos_x << " " << pos_y << "\n";
	if (pos_x + 1 == pos_y)
	{
		// std::cout << "A\n";
		curr_delta = -data->time_cost[prev_x][cl_x] + data->time_cost[prev_x][cl_y] -
			data->time_cost[cl_y][cl_w] + data->time_cost[cl_x][cl_w];
	}
	// else if (pos_y + 1 == pos_x)
	// {
	// 	// std::cout << "A\n";
	// 	curr_delta = -data->time_cost[prev_y][cl_y] + data->time_cost[prev_y][cl_x] -
	// 		data->time_cost[cl_x][cl_z] + data->time_cost[cl_y][cl_z];
	// }
	else
		curr_delta = -data->time_cost[prev_x][cl_x] - data->time_cost[cl_x][cl_z] + data->time_cost[prev_x][cl_y] + data->time_cost[cl_y][cl_z] +
			-data->time_cost[prev_y][cl_y] - data->time_cost[cl_y][cl_w] + data->time_cost[prev_y][cl_x] + data->time_cost[cl_x][cl_w];

#ifdef TW
	if (route_x->is_time_feasible && curr_delta >= EPSILON)
		return;

	if (pos_x < pos_y)
	{
		ConcatRouteXSwap11IntraFront();

#ifdef ROBUST_TW
		curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation);
#else
		curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
#endif
	}
	else
	{
		ConcatRouteXSwap11IntraBack();
#ifdef ROBUST_TW
		curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation);
#else
		curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
#endif
	}
#endif

#ifdef USEMMD
	if (curr_delta + EPSILON < intra_smd->cost)
	{
		intra_smd->status = IMPROVED;
		intra_smd->cost = curr_delta;
		intra_smd->r1_pos1 = pos_x;
		intra_smd->r2_pos1 = pos_y;
	}

#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
	}
#endif
}

inline void LocalSearch::ExecSwap11()
{
	std::swap(curr_solution->route[best_route_x->cour][best_pos_x], curr_solution->route[best_route_y->cour][best_pos_y]);
	if (best_route_x != best_route_y)
	{
		best_route_x->Update();
		best_route_y->Update();
	}
	else
		best_route_x->Update();
	UpdateClientRank();
}

inline void LocalSearch::SearchSwap11()
{
	//     improved = false;
	//     best_delta = 0;
	//     for (cl_x = 1; cl_x <= data->nb_clients; cl_x++)
	//     {
	//         for (int cy = 0; cy < data->correlated_vertices[cl_x].size(); cy++)
	//         {
	//             // prev_y = data->correlated_vertices[cl_x][cy];
	//             // routeY = clientRank[prev_y].first;
	//             // pos_y = clientRank[prev_y].second + 1;
	//             // cl_y = curr_solution->route[ry_idx][pos_y];
	// 
	//             // if (pos_y == route_y->length - 1)
	//             //     continue;
	// 
	//             // pos_x = clientRank[cl_x].second;
	//             // pos_z = pos_x + 1;
	//             // pos_w = pos_y + 1;
	//             // routeX = clientRank[cl_x].first;
	//             // cl_z = curr_solution->route[rx_idx][pos_z];
	//             // cl_w = curr_solution->route[ry_idx][pos_w];
	//             // prev_x = curr_solution->route[rx_idx][pos_x - 1];
	// 
	//             // if (cl_x == cl_y)
	//             //     continue;
	// 
	//             // // if(pos_x + 1 == pos_y)
	//             // //             std::cout << "A\n";
	// 
	//             // // SetLocalVariables();
	// 
	//             cl_y = data->correlated_vertices[cl_x][cy];
	//             SetLocalVariables();
	//             EvalSwap11();
	//         }
	//     }
	//     if (improved)
	//     {
	//         // std::cout << best_delta << "\n";
	//         // std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";
	//         ExecSwap11();
	//     }
}

void LocalSearch::SearchSwap11Complete()
{
	stats->inter_neighb_explore[1].SetStart();
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		route_x = ads->routeADS[rx_idx];
		if (route_x->length <= 2)
			continue;

		for (ry_idx = rx_idx + 1; ry_idx < ads->routeADS.size(); ry_idx++)
		{
			route_y = ads->routeADS[ry_idx];
			rx_ry_feasible = route_x->is_feasible && route_y->is_feasible;

			if (route_y->length <= 2
					|| (rx_ry_feasible && (route_y->load.back() - route_y->max_demand + route_x->min_demand > data->capacity + EPSILON))
					|| (!curr_solution->GetInterNeighbStatus(rx_idx, ry_idx, 1)))
				continue;

			rx_ry_improved  = false;

			stats->inter_mmd_search[1].SetStart();
#ifdef USEMMD
			inter_smd 		= ads->getSMD(rx_idx, ry_idx, 1);	
			inter_smd_rev = NULL;
#endif
			stats->inter_mmd_search[1].SetEnd();

			start_search = 0;
			end_search  = route_y->length - 2;

			if (!use_mmd || inter_smd->status == UNVISITED)
			{
#ifdef USEMMD
				inter_smd->status = UNIMPROVED;
				inter_smd_rev = ads->getSMD(ry_idx, rx_idx, 1);	
				inter_smd_rev->status = UNIMPROVED;
#endif

				for (pos_x = 1; pos_x < route_x->length - 1; pos_x++)
				{
					cl_x = curr_solution->route[rx_idx][pos_x];
					if (rx_ry_feasible && (route_y->load.back() - route_y->max_demand + data->client[cl_x].demand > data->capacity))
						continue;
					prev_x = curr_solution->route[rx_idx][pos_x - 1];
					pos_z = pos_x + 1;
					cl_z = curr_solution->route[rx_idx][pos_z];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_x];
						end_search = route_y->end_search[cl_x];
						// std::cout << "\n";
						// std::cout << start_search << " " << 0 << "\n";
						// std::cout << end_search << " " << route_y->length - 2 << "\n";
						// std::cout << "\n";
					}
#endif

					for (pos_y = start_search + 1; pos_y <= end_search; pos_y++)
					{
						pos_w = pos_y + 1;
						cl_y = curr_solution->route[ry_idx][pos_y];
						cl_w = curr_solution->route[ry_idx][pos_w];
						prev_y = curr_solution->route[ry_idx][pos_y - (pos_y > 0)];

						EvalSwap11();
						rx_ry_improved |= (curr_delta + EPSILON < 0) ;
					}
				}
			}	else 
			{
				stats->inter_mmd_effective[1]++;
			}

#ifdef USEMMD		
			rx_ry_improved |= (inter_smd->cost + EPSILON < 0);

			if (inter_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = inter_smd->cost;
				best_pos_x = inter_smd->r1_pos1;
				best_pos_y = inter_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
			}

			if (inter_smd_rev)
				*inter_smd_rev = inter_smd->Reverse();
#endif

			if (!rx_ry_improved)	
				curr_solution->SetInterNeighbStatus(rx_idx, ry_idx, 1, false);
		}
	}

	if (improved)
	{
		stats->neighb_improv[2]++;
		// std::cout << *curr_solution << "\n";
		// std::cout << best_delta << "\n";
		// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";

		curr_solution->SetInterNeighbStatusAllTrue(best_route_x->cour, best_route_y->cour);
		ExecSwap11();
	}
	stats->inter_neighb_explore[1].SetEnd();
}

/** ********************* SWAP (2,2) ***************************** */

inline void LocalSearch::ComputeRobustDemandRouteXSwap22()
{
	routeXRobustDemand = route_x->load.back() - data->client[cl_x].demand -
		data->client[cl_z].demand + data->client[cl_y].demand + data->client[cl_w].demand;
	std::vector<double> sorted_deviations(route_x->length - 2);

	int i;
	for (i = 1; i <= pos_x - 1; i++)
	{
		sorted_deviations[i - 1] = -data->client[curr_solution->route[rx_idx][i]].demand_deviation;
	}
	sorted_deviations[i++ - 1] = -data->client[cl_y].demand_deviation;
	sorted_deviations[i++ - 1] = -data->client[cl_w].demand_deviation;
	for (; i <= route_x->length - 2; i++)
	{
		sorted_deviations[i - 1] = -data->client[curr_solution->route[rx_idx][i]].demand_deviation;
	}

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}
	routeXRobustCapViolation = std::max(0.0, routeXRobustDemand - data->capacity);
}

inline void LocalSearch::ComputeRobustDemandRouteYSwap22()
{
	routeYRobustDemand = route_y->load.back() + data->client[cl_x].demand +
		data->client[cl_z].demand - data->client[cl_y].demand - data->client[cl_w].demand;
	std::vector<double> sorted_deviations(route_y->length - 2);

	int i;
	for (i = 1; i <= pos_y - 1; i++)
	{
		sorted_deviations[i - 1] = -data->client[curr_solution->route[ry_idx][i]].demand_deviation;
	}
	sorted_deviations[i++ - 1] = -data->client[cl_x].demand_deviation;
	sorted_deviations[i++ - 1] = -data->client[cl_z].demand_deviation;
	for (; i <= route_y->length - 2; i++)
	{
		sorted_deviations[i - 1] = -data->client[curr_solution->route[ry_idx][i]].demand_deviation;
	}

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}

	routeYRobustCapViolation = std::max(0.0, routeYRobustDemand - data->capacity);
}

inline void LocalSearch::ConcatRouteXSwap22Inter()
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	pathX.push_back(cl_w);
	for (int i = pos_x + 2; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
	routeXRobustTWViolationRev = 0;
	std::swap(pathX[pos_x+1], pathX[pos_x]);
	ComputePathEarliestStart(pathX, routeXRobustTWViolationRev);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_y->subseq[pos_y][pos_y + 1]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 2][route_x->length - 1]);
	subseq_x_rev = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_y->subseq[pos_y + 1][pos_y]);
	subseq_x_rev = Subsequence::Concatenate(data, subseq_x_rev, route_x->subseq[pos_x + 2][route_x->length - 1]);
}

inline void LocalSearch::ConcatRouteYSwap22Inter()
{
#ifdef ROBUST_TW
	routeYRobustTWViolation = 0;
	pathY.clear();
	for (int i = 0; i <= pos_y - 1; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	pathY.push_back(cl_x);
	pathY.push_back(cl_z);
	for (int i = pos_y + 2; i < route_y->length; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	ComputePathEarliestStart(pathY, routeYRobustTWViolation);
	routeYRobustTWViolationRev = 0;
	std::swap(pathY[pos_y+1], pathY[pos_y]);
	ComputePathEarliestStart(pathY, routeYRobustTWViolationRev);
#endif

	subseq_y = Subsequence::Concatenate(data, route_y->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x + 1]);
	subseq_y = Subsequence::Concatenate(data, subseq_y, route_y->subseq[pos_y + 2][route_y->length - 1]);
	subseq_y_rev = Subsequence::Concatenate(data, route_y->subseq[0][pos_y - 1], route_x->subseq[pos_x + 1][pos_x]);
	subseq_y_rev = Subsequence::Concatenate(data, subseq_y_rev, route_y->subseq[pos_y + 2][route_y->length - 1]);
}

inline void LocalSearch::ConcatRouteXSwap22IntraFront() // when pos_x < pos_y
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	pathX.push_back(cl_w);
	for (int i = pos_x + 2; i <= pos_y - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_x);
	pathX.push_back(cl_z);
	for (int i = pos_y + 2; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif
	if (pos_z + 1 == pos_y)
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_y][pos_y + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x][pos_x + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 2][route_x->length - 1]);
	}
	else
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_y][pos_y + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 2][pos_y - 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x][pos_x + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 2][route_x->length - 1]);
	}
}

inline void LocalSearch::ConcatRouteXSwap22IntraBack() // when pos_x > pos_y
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_y - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_x);
	pathX.push_back(cl_z);
	for (int i = pos_y + 2; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	pathX.push_back(cl_w);
	for (int i = pos_x + 2; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	if (pos_w + 1 == pos_x)
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y][pos_y + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 2][route_x->length - 1]);
	}
	else
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 2][pos_x - 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y][pos_y + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 2][route_x->length - 1]);
	}
}

inline void LocalSearch::EvalSwap22()
{
	double temp = temp_delta
								- data->time_cost[prev_y][cl_y] 
								- data->time_cost[cl_w][next_w];

	curr_delta1 = data->time_cost[prev_x][cl_y] + data->time_cost[cl_w][next_z] +
		+data->time_cost[prev_y][cl_x] + data->time_cost[cl_z][next_w] + temp;

	curr_delta2 = data->time_cost[prev_x][cl_y] + data->time_cost[cl_w][next_z] +
		+data->time_cost[prev_y][cl_z] + data->time_cost[cl_x][next_w] + temp;

	curr_delta3 = data->time_cost[prev_x][cl_w] + data->time_cost[cl_y][next_z] +
		+data->time_cost[prev_y][cl_x] + data->time_cost[cl_z][next_w] + temp;

	curr_delta4 = data->time_cost[prev_x][cl_w] + data->time_cost[cl_y][next_z] +
		+data->time_cost[prev_y][cl_z] + data->time_cost[cl_x][next_w] + temp;
		
	curr_delta=0;

	if (route_x->is_feasible && route_y->is_feasible && curr_delta1 >= EPSILON && curr_delta2 >= EPSILON && curr_delta3 >= EPSILON && curr_delta4 >= EPSILON)
		return;

#ifdef TW
	ConcatRouteXSwap22Inter();
	ConcatRouteYSwap22Inter();

#ifdef ROBUST_TW
	curr_delta1 += 1 * data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		1 * data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolation);
	curr_delta2 += 1 * data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		1 * data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolationRev);
	curr_delta3 += 1 * data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolationRev) +
		1 * data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolation);
	curr_delta4 += 1 * data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolationRev) +
		1 * data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolationRev);
#else
	curr_delta1 += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp) +
		data->penalty_tw * (-route_y->tw_violation + subseq_y.timewarp);
	curr_delta2 += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp) +
		data->penalty_tw * (-route_y->tw_violation + subseq_y_rev.timewarp);
	curr_delta3 += data->penalty_tw * (-route_x->tw_violation + subseq_x_rev.timewarp) +
		data->penalty_tw * (-route_y->tw_violation + subseq_y.timewarp);
	curr_delta4 += data->penalty_tw * (-route_x->tw_violation + subseq_x_rev.timewarp) +
		data->penalty_tw * (-route_y->tw_violation + subseq_y_rev.timewarp);
#endif
#endif

#ifdef ROBUST_DEMAND
	ComputeRobustDemandRouteXSwap22();
	ComputeRobustDemandRouteYSwap22();
	temp = data->penalty_load * (routeXRobustCapViolation + routeYRobustCapViolation -
			route_x->robust_cap_violation - route_y->robust_cap_violation);
	curr_delta1 += temp;
	curr_delta2 += temp;
	curr_delta3 += temp;
	curr_delta4 += temp;
#else
	double tempLoadDelta = -data->client[cl_x].demand - data->client[cl_z].demand + data->client[cl_y].demand + data->client[cl_w].demand;
	temp = data->penalty_load * (-route_x->cap_violation + ComputeCapViolation(route_x->load.back() + tempLoadDelta)) +
		data->penalty_load * (-route_y->cap_violation + ComputeCapViolation(route_y->load.back() - tempLoadDelta));

	if (rx_ry_feasible && temp > 0) {
		return;
	}

	//curr_delta1 += temp;
	//curr_delta2 += temp;
	//curr_delta3 += temp;
	//curr_delta4 += temp;
#endif

	curr_delta = curr_delta1;
	move_type = 0;
	if (curr_delta2 < curr_delta)
	{
		curr_delta = curr_delta2;
		move_type = 2;
	}
	if (curr_delta3 < curr_delta)
	{
		curr_delta = curr_delta3;
		move_type = 3;
	}
	if (curr_delta4 < curr_delta)
	{
		curr_delta = curr_delta4;
		move_type = 1;
	}

	curr_delta += temp;

#ifdef USEMMD
	if (curr_delta + EPSILON < inter_smd->cost)
	{
		inter_smd->status = IMPROVED;
		inter_smd->cost = curr_delta;
		inter_smd->r1_pos1 = pos_x;
		inter_smd->r2_pos1 = pos_y;
		inter_smd->move_type = move_type;
	}
#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
		bestMoveType = move_type;
	}
#endif
}

// inline void LocalSearch::EvalExchange22()
// {
// 	if (pos_z + 1 == pos_y)
// 	{
// 		curr_delta = (-data->time_cost[prev_x][cl_x] + data->time_cost[prev_x][cl_y] -
// 				data->time_cost[cl_w][next_w] + data->time_cost[cl_z][next_w]) -
// 			data->time_cost[cl_z][cl_y] + data->time_cost[cl_w][cl_x];
// 	}
// 	else if (pos_w + 1 == pos_x)
// 	{
// 		// return;
// 		curr_delta = (-data->time_cost[prev_y][cl_y] + data->time_cost[prev_y][cl_x] +
// 				-data->time_cost[cl_z][next_z] + data->time_cost[cl_w][next_z]) -
// 			data->time_cost[cl_w][cl_x] + data->time_cost[cl_y][cl_z];
// 	}
// 	else
// 		curr_delta = (-data->time_cost[prev_x][cl_x] - data->time_cost[cl_z][next_z] + data->time_cost[prev_x][cl_y] + data->time_cost[cl_w][next_z] +
// 				-data->time_cost[prev_y][cl_y] - data->time_cost[cl_w][next_w] + data->time_cost[prev_y][cl_x] + data->time_cost[cl_z][next_w]);
// 
// #ifdef TW
// 	if (route_x->is_time_feasible && curr_delta >= EPSILON)
// 		return;
// 
// 	if (pos_x < pos_y)
// 	{
// 		ConcatRouteXSwap22IntraFront();
// #ifdef ROBUST_TW
// 		curr_delta += 1 * data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation);
// #else
// 		curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
// #endif
// 	}
// 	else
// 	{
// 		ConcatRouteXSwap22IntraBack();
// #ifdef ROBUST_TW
// 		curr_delta += 1 * data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation);
// #else
// 		curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
// #endif
// 	}
// #endif
// 
// 	if (curr_delta + EPSILON < best_delta)
// 	{
// 		improved = true;
// 		best_delta = curr_delta;
// 		best_pos_x = pos_x;
// 		best_pos_y = pos_y;
// 		best_route_x = route_x;
// 		best_route_y = route_y;
// 	}
// }

inline void LocalSearch::ExecSwap22()
{
	std::swap(curr_solution->route[best_route_x->cour][best_pos_x], curr_solution->route[best_route_y->cour][best_pos_y]);
	std::swap(curr_solution->route[best_route_x->cour][best_pos_x + 1], curr_solution->route[best_route_y->cour][best_pos_y + 1]);
	if (bestMoveType == 2 || bestMoveType == 1)
	{
		std::swap(curr_solution->route[best_route_y->cour][best_pos_y], curr_solution->route[best_route_y->cour][best_pos_y + 1]);
	}
	if (bestMoveType == 3 || bestMoveType == 1)
	{
		std::swap(curr_solution->route[best_route_x->cour][best_pos_x], curr_solution->route[best_route_x->cour][best_pos_x + 1]);
	}

	if (best_route_x != best_route_y)
	{
		best_route_x->Update();
		best_route_y->Update();
	}
	else
		best_route_x->Update();

	UpdateClientRank();
}

inline void LocalSearch::SearchSwap22()
{
	//     improved = false;
	//     best_delta = 0;
	//     for (cl_x = 1; cl_x <= data->nb_clients; cl_x++)
	//     {
	//         for (int cy = 0; cy < data->correlated_vertices[cl_x].size(); cy++)
	//         {
	//             prev_y = data->correlated_vertices[cl_x][cy];
	//             routeY = clientRank[prev_y].first;
	//             pos_y = clientRank[prev_y].second + 1;
	//             cl_y = curr_solution->route[ry_idx][pos_y];
	// 
	//             if (cl_y == cl_x)
	//                 continue;
	// 
	//             if (pos_y >= route_y->length - 3 || route_y->length < 4)
	//                 continue;
	// 
	//             pos_x = clientRank[cl_x].second;
	//             pos_z = pos_x + 1;
	//             pos_w = pos_y + 1;
	//             routeX = clientRank[cl_x].first;
	//             if (pos_x >= route_x->length - 3 || route_x->length < 4)
	//                 continue;
	// 
	//             cl_z = curr_solution->route[rx_idx][pos_z];
	//             cl_w = curr_solution->route[ry_idx][pos_w];
	//             prev_x = curr_solution->route[rx_idx][pos_x - 1];
	//             next_w = curr_solution->route[ry_idx][pos_w + 1];
	//             next_z = curr_solution->route[rx_idx][pos_z + 1];
	// 
	//             // if (routeX != routeY)
	//             //     continue;
	// 
	//             // SetLocalVariables();
	//             EvalSwap22();
	//         }
	//     }
	//     if (improved)
	//     {
	//         // std::cout << best_delta << "\n";
	//         // std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";
	//         ExecSwap22();
	//     }
}

void LocalSearch::SearchSwap22Complete()
{
	stats->inter_neighb_explore[2].SetStart();
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		route_x = ads->routeADS[rx_idx];
		if (route_x->length <= 3)
			continue;

		for (ry_idx = rx_idx + 1; ry_idx < ads->routeADS.size(); ry_idx++)
		{
			route_y = ads->routeADS[ry_idx];
			rx_ry_feasible = route_x->is_feasible && route_y->is_feasible;
			rx_ry_improved = false;

			if (route_y->length <= 3
					|| (rx_ry_feasible && (route_y->load.back() - route_y->max_pair_demand + route_x->min_pair_demand > data->capacity + EPSILON))
				  || (!curr_solution->GetInterNeighbStatus(rx_idx, ry_idx, 2)))
				continue;

			stats->inter_mmd_search[2].SetStart();
#ifdef USEMMD
			inter_smd = ads->getSMD(rx_idx, ry_idx, 2);	
			inter_smd_rev = NULL;
#endif
			stats->inter_mmd_search[2].SetEnd();

			start_search = 0;
			end_search  = route_y->length - 2;

			if (!use_mmd || inter_smd->status == UNVISITED)
			{
#ifdef USEMMD
				inter_smd->status = UNIMPROVED;
				inter_smd_rev = ads->getSMD(ry_idx, rx_idx, 2);	
				inter_smd_rev->status = UNIMPROVED;
#endif

				for (pos_x = 1; pos_x < route_x->length - 2; pos_x++)
				{
					pos_z = pos_x + 1;
					cl_x = curr_solution->route[rx_idx][pos_x];
					cl_z = curr_solution->route[rx_idx][pos_z];
					if (rx_ry_feasible && (route_y->load.back() - route_y->max_pair_demand
																											+ data->client[cl_x].demand 
																											+ data->client[cl_z].demand > data->capacity))
						continue;

					prev_x = curr_solution->route[rx_idx][pos_x - 1];
					next_z = curr_solution->route[rx_idx][pos_z + 1];
					temp_delta = -data->time_cost[prev_x][cl_x] - data->time_cost[cl_z][next_z];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_x];
						end_search = route_y->end_search[cl_x];
						// std::cout << "\n";
						// std::cout << start_search << " " << 0 << "\n";
						// std::cout << end_search << " " << route_y->length - 2 << "\n";
						// std::cout << "\n";
					}
#endif

					for (pos_y = start_search + 1; pos_y <= end_search - 1; pos_y++)
					{
						pos_w  = pos_y + 1;
						cl_y   = curr_solution->route[ry_idx][pos_y];
						cl_w   = curr_solution->route[ry_idx][pos_w];
						prev_y = curr_solution->route[ry_idx][pos_y - 1];
						next_w = curr_solution->route[ry_idx][pos_w + 1];

						EvalSwap22();
						rx_ry_improved |= (curr_delta + EPSILON < 0);
					}
				}
			}
			else
			{
				stats->inter_mmd_effective[2]++;
			}
#ifdef USEMMD		
			rx_ry_improved |= (inter_smd->cost + EPSILON < 0);

			if (inter_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = inter_smd->cost;
				best_pos_x = inter_smd->r1_pos1;
				best_pos_y = inter_smd->r2_pos1;
				bestMoveType = inter_smd->move_type;
				best_route_x = route_x;
				best_route_y = route_y;
			}

			if (inter_smd_rev)
				*inter_smd_rev = inter_smd->Reverse();
#endif

			 if (!rx_ry_improved)
				 curr_solution->SetInterNeighbStatus(rx_idx, ry_idx, 2, false);
		}
	}

	if (improved)
	{
		stats->neighb_improv[4]++;
		// std::cout << *curr_solution << "\n";
		// std::cout << best_delta << "\n";
		// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";

		curr_solution->SetInterNeighbStatusAllTrue(best_route_x->cour, best_route_y->cour);
		ExecSwap22();
	}
	stats->inter_neighb_explore[2].SetEnd();
}

/** ********************* CROSS ***************************** */

inline void LocalSearch::ComputeRobustDemandRouteXCross()
{
	routeXRobustDemand = route_x->load[pos_x] + route_y->load.back() - route_y->load[pos_y];
	std::vector<double> sorted_deviations;

	int i;
	for (i = 1; i <= pos_x; i++)
	{
		sorted_deviations.push_back(-data->client[curr_solution->route[rx_idx][i]].demand_deviation);
	}

	for (int j = pos_w; j <= route_y->length - 2; i++, j++)
	{
		sorted_deviations.push_back(-data->client[curr_solution->route[ry_idx][j]].demand_deviation);
	}

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}
	routeXRobustCapViolation = std::max(0.0, routeXRobustDemand - data->capacity);
}

inline void LocalSearch::ComputeRobustDemandRouteYCross()
{
	double routeYRobustDemand = route_y->load[pos_y] + route_x->load.back() - route_x->load[pos_x];
	std::vector<double> sorted_deviations;

	int i;
	for (i = 1; i <= pos_y; i++)
	{
		sorted_deviations.push_back(-data->client[curr_solution->route[ry_idx][i]].demand_deviation);
	}
	for (int j = pos_z; j <= route_x->length - 2; i++, j++)
	{
		sorted_deviations.push_back(-data->client[curr_solution->route[rx_idx][j]].demand_deviation);
	}

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}

	routeYRobustCapViolation = std::max(0.0, routeYRobustDemand - data->capacity);
}

inline void LocalSearch::ConcatRouteXCross()
{
#ifdef ROBUST_TW
	pathX.clear();
	routeXRobustTWViolation = 0;
	for (int i = 0; i <= pos_x; i++)
		pathX.push_back(curr_solution->route[rx_idx][i]);

	for (int i = pos_x + 1, j = pos_y + 1; j < route_y->length; j++, i++)
		pathX.push_back(curr_solution->route[ry_idx][j]);

	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x], route_y->subseq[pos_w][route_y->length - 1]);
}

inline void LocalSearch::ConcatRouteYCross()
{
#ifdef ROBUST_TW
	pathY.clear();
	routeYRobustTWViolation = 0;
	for (int i = 0; i <= pos_y; i++)
		pathY.push_back(curr_solution->route[ry_idx][i]);

	for (int i = pos_y + 1, j = pos_x + 1; j < route_x->length; j++, i++)
		pathY.push_back(curr_solution->route[rx_idx][j]);

	ComputePathEarliestStart(pathY, routeYRobustTWViolation);
#endif

	subseq_y = Subsequence::Concatenate(data, route_y->subseq[0][pos_y], route_x->subseq[pos_z][route_x->length - 1]);
}

inline void LocalSearch::EvalCross()
{
	curr_delta = -data->time_cost[cl_x][cl_z] - data->time_cost[cl_y][cl_w] +
		data->time_cost[cl_x][cl_w] + data->time_cost[cl_y][cl_z];
	// curr_delta = -route_x->distance + subseq_x.distance +
	//             (-route_y->distance + subseq_y.distance);

	if (route_x->is_feasible && route_y->is_feasible && curr_delta >= EPSILON)
		return;

#ifdef TW
	ConcatRouteXCross();
	ConcatRouteYCross();

#ifdef ROBUST_TW

	curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolation);

#else
	curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp) +
		data->penalty_tw * (-route_y->tw_violation + subseq_y.timewarp);

#endif
#endif

#ifdef ROBUST_DEMAND
	ComputeRobustDemandRouteXCross();
	ComputeRobustDemandRouteYCross();
	curr_delta += data->penalty_load * (routeXRobustCapViolation + routeYRobustCapViolation -
			route_x->robust_cap_violation - route_y->robust_cap_violation);
#else
	curr_delta += data->penalty_load * (-route_x->cap_violation + ComputeCapViolation(route_x->load[pos_x] + route_y->load.back() - route_y->load[pos_y])) +
		data->penalty_load * (-route_y->cap_violation + ComputeCapViolation(route_y->load[pos_y] + route_x->load.back() - route_x->load[pos_x]));
#endif

#ifdef USEMMD
	if (curr_delta + EPSILON < inter_smd->cost)
	{
		inter_smd->status = IMPROVED;
		inter_smd->cost = curr_delta;
		inter_smd->r1_pos1 = pos_x;
		inter_smd->r2_pos1 = pos_y;
	}
#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
	}
#endif
}

inline void LocalSearch::ExecCross()
{
	std::vector<int> v1(curr_solution->route[best_route_x->cour].begin() + best_pos_x + 1, curr_solution->route[best_route_x->cour].end()), v2(curr_solution->route[best_route_y->cour].begin() + best_pos_y + 1, curr_solution->route[best_route_y->cour].end());
	curr_solution->route[best_route_x->cour].erase(curr_solution->route[best_route_x->cour].begin() + best_pos_x + 1, curr_solution->route[best_route_x->cour].end());
	curr_solution->route[best_route_y->cour].erase(curr_solution->route[best_route_y->cour].begin() + best_pos_y + 1, curr_solution->route[best_route_y->cour].end());
	curr_solution->route[best_route_x->cour].insert(curr_solution->route[best_route_x->cour].end(), v2.begin(), v2.end());
	curr_solution->route[best_route_y->cour].insert(curr_solution->route[best_route_y->cour].end(), v1.begin(), v1.end());

	best_route_x->Update();
	best_route_y->Update();

	UpdateClientRank();
}

inline void LocalSearch::SearchCross()
{
	//     improved = false;
	//     best_delta = 0;
	//     for (cl_x = 1; cl_x <= data->nb_clients; cl_x++)
	//     {
	//         for (int cw = 0; cw < data->correlated_vertices[cl_x].size(); cw++)
	//         {
	//             cl_w = data->correlated_vertices[cl_x][cw];
	//             pos_w = clientRank[cl_w].second;
	//             routeY = clientRank[cl_w].first;
	//             pos_y = pos_w - 1;
	//             // prev_y = curr_solution->route[ry_idx][pos_y - 1];
	//             cl_y = curr_solution->route[ry_idx][pos_y];
	// 
	//             pos_x = clientRank[cl_x].second;
	//             pos_z = pos_x + 1;
	//             pos_w = pos_y + 1;
	//             routeX = clientRank[cl_x].first;
	//             cl_z = curr_solution->route[rx_idx][pos_z];
	//             prev_x = curr_solution->route[rx_idx][pos_x - 1];
	// 
	//             // SetLocalVariables();
	//             EvalCross();
	//         }
	//     }
	//     if (improved)
	//     {
	//         // std::cout << "\n";
	//         // std::cout << *curr_solution << "\n";
	//         // std::cout << best_delta << "\n";
	//         // std::cout << best_pos_x << '\n';
	//         // std::cout << best_pos_y << '\n';
	//         // std::cout << best_route_x << '\n';
	//         // std::cout << best_route_y << '\n';
	//         // std::cout << "\n";
	//         // std::cout << best_delta << "\n";
	//         ExecCross();
	//     }
}

void LocalSearch::SearchCrossComplete()
{
	stats->inter_neighb_explore[3].SetStart();
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		route_x = ads->routeADS[rx_idx];
		if (route_x->length <= 2)
			continue;

		// for (ry_idx = (rx_idx + 1); ry_idx < ads->routeADS.size(); ry_idx++)
		for (ry_idx = 0; ry_idx < ads->routeADS.size(); ry_idx++)
		{
			route_y = ads->routeADS[ry_idx];
			if (rx_idx == ry_idx || route_y->length <= 2 
				 || (!curr_solution->GetInterNeighbStatus(rx_idx, ry_idx, 3)))
				continue;

			stats->inter_mmd_search[3].SetStart();
#ifdef USEMMD
			inter_smd = ads->getSMD(rx_idx, ry_idx, 3);	
			inter_smd_rev = NULL;
#endif
			stats->inter_mmd_search[3].SetEnd();

			start_search = 0;
			end_search  = route_y->length - 2;
			rx_ry_feasible = route_x->is_feasible && route_y->is_feasible;
			rx_ry_improved = false;

			if (!use_mmd || inter_smd->status == UNVISITED)
			{
#ifdef USEMMD
				inter_smd->status = UNIMPROVED;
				inter_smd_rev = ads->getSMD(ry_idx, rx_idx, 3);	
				inter_smd_rev->status = UNIMPROVED;
#endif

				for (pos_x = 1; pos_x < route_x->length - 1; pos_x++)
				{
					// if (rx_ry_feasible && (route_x->load[pos_x] + 
					// 			data->client[curr_solution->route[ry_idx].demand[curr_solution->route[ry_idx].size() - 2]] 
					// 			> data->capacity))
					// 	continue;

					cl_x = curr_solution->route[rx_idx][pos_x];
					pos_z = pos_x + 1;
					cl_z = curr_solution->route[rx_idx][pos_z];
					// prev_x = curr_solution->route[rx_idx][pos_x - 1];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_z];
						end_search = route_y->end_search[cl_z];
						// std::cout << "\n";
						// std::cout << start_search << " " << 0 << "\n";
						// std::cout << end_search << " " << route_y->length - 2 << "\n";
						// std::cout << "\n";
					}
#endif

					for (pos_y = end_search; pos_y >= start_search; pos_y--)
					{
						pos_w = pos_y + 1;
						cl_y = curr_solution->route[ry_idx][pos_y];
						cl_w = curr_solution->route[ry_idx][pos_w];

						if (rx_ry_feasible && 
								(route_x->load[pos_x] - route_y->load[pos_y] + route_y->load.back() > data->capacity))
							break;

						EvalCross();
						rx_ry_improved |= (curr_delta + EPSILON < 0);
					}
				}
			}
			else
			{
				stats->inter_mmd_effective[3]++;
			}
#ifdef USEMMD		
			rx_ry_improved |= (inter_smd->cost + EPSILON < 0);

			if (inter_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = inter_smd->cost;
				best_pos_x = inter_smd->r1_pos1;
				best_pos_y = inter_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
			}

			if (inter_smd_rev)
				*inter_smd_rev = inter_smd->Reverse();
#endif

			 if (!rx_ry_improved)
			 	curr_solution->SetInterNeighbStatus(rx_idx, ry_idx, 3, false);
		}
	}

	if (improved)
	{
		stats->neighb_improv[5]++;
		// std::cout << *curr_solution << "\n";
		// std::cout << best_delta << "\n";
		// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x->cour << " " << best_route_y->cour << "\n";

		curr_solution->SetInterNeighbStatusAllTrue(best_route_x->cour, best_route_y->cour);
		ExecCross();
	}
	stats->inter_neighb_explore[3].SetEnd();
}

/** ********************* RELOCATE 2 ***************************** */

inline void LocalSearch::ComputeRobustDemandRouteXRelocate2()
{
	routeXRobustDemand = route_x->load.back() - data->client[cl_x].demand - data->client[cl_z].demand;
	std::vector<double> sorted_deviations(route_x->length - 4);

	for (int i = 1, j = 0; i <= route_x->length - 2; i++)
	{
		if (i == pos_x || i == pos_z)
			continue;
		sorted_deviations[j] = -data->client[curr_solution->route[rx_idx][i]].demand_deviation;
		j++;
	}

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}

	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}

	routeXRobustCapViolation = std::max(0.0, routeXRobustDemand - data->capacity);
}

inline void LocalSearch::ComputeRobustDemandRouteYRelocate2()
{
	double routeYRobustDemand = route_y->load.back() + data->client[cl_x].demand + data->client[cl_z].demand;
	std::vector<double> sorted_deviations(route_y->length);

	int i;
	for (i = 1; i <= route_y->length - 2; i++)
		sorted_deviations[i - 1] = -data->client[curr_solution->route[ry_idx][i]].demand_deviation;
	sorted_deviations.back() = -data->client[cl_x].demand_deviation;
	sorted_deviations[sorted_deviations.size() - 2] = -data->client[cl_z].demand_deviation;

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}

	routeYRobustCapViolation = std::max(0.0, routeYRobustDemand - data->capacity);
}

inline void LocalSearch::ConcatRouteXRelocate2Inter()
{
#ifdef ROBUST_TW
	// int l = 0;
	// ComputeNodeEarliestStartRouteX(curr_solution->route[rx_idx][0], 0, l++);

	// for (int j = 1; j <= pos_x - 1; j++)
	//     ComputeNodeEarliestStartRouteX(curr_solution->route[rx_idx][j], curr_solution->route[rx_idx][j - 1], l++);

	// ComputeNodeEarliestStartRouteX(curr_solution->route[rx_idx][pos_x + 1], curr_solution->route[rx_idx][pos_x - 1], l++);

	// for (int j = pos_x + 2; j <= route_x->length - 1; j++)
	//     ComputeNodeEarliestStartRouteX(curr_solution->route[rx_idx][j], curr_solution->route[rx_idx][j - 1], l++);

	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_x + 2; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_x + 2][route_x->length - 1]);
}

inline void LocalSearch::ConcatRouteYRelocate2Inter()
{
#ifdef ROBUST_TW
	// int l = 0;
	// ComputeNodeEarliestStartRouteY(curr_solution->route[ry_idx][0], 0, l++);
	// for (int j = 1; j <= pos_y; j++)
	//     ComputeNodeEarliestStartRouteY(curr_solution->route[ry_idx][j], curr_solution->route[ry_idx][j - 1], l++);

	// ComputeNodeEarliestStartRouteY(curr_solution->route[rx_idx][pos_x], curr_solution->route[ry_idx][pos_y], l++);
	// ComputeNodeEarliestStartRouteY(curr_solution->route[ry_idx][pos_y + 1], curr_solution->route[rx_idx][pos_x], l++);

	// for (int j = pos_y + 2; j <= route_y->length - 1; j++)
	//     ComputeNodeEarliestStartRouteY(curr_solution->route[ry_idx][j], curr_solution->route[ry_idx][j - 1], l++);
	routeYRobustTWViolation = 0;
	pathY.clear();
	for (int i = 0; i <= pos_y; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	pathY.push_back(curr_solution->route[rx_idx][pos_x]);
	pathY.push_back(curr_solution->route[rx_idx][pos_x + 1]);
	for (int i = pos_y + 1; i < route_y->length; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	ComputePathEarliestStart(pathY, routeYRobustTWViolation);
	routeYRobustTWViolationRev = 0;
	std::swap(pathY[pos_y+2], pathY[pos_y+1]);
	ComputePathEarliestStart(pathY, routeYRobustTWViolationRev);
#endif

	subseq_y = Subsequence::Concatenate(data, route_y->subseq[0][pos_y], route_x->subseq[pos_x][pos_x + 1]);
	subseq_y = Subsequence::Concatenate(data, subseq_y, route_y->subseq[pos_y + 1][route_y->length - 1]);
	subseq_y_rev = Subsequence::Concatenate(data, route_y->subseq[0][pos_y], route_x->subseq[pos_x + 1][pos_x]);
	subseq_y_rev = Subsequence::Concatenate(data, subseq_y_rev, route_y->subseq[pos_y + 1][route_y->length - 1]);
}

inline void LocalSearch::ConcatRouteXRelocate2IntraFront()
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_z + 1; i <= pos_y; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	//pathX.push_back(cl_x);
	//pathX.push_back(cl_z);
	for (int i = pos_x; i <= pos_z; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_y + 1; i <= route_x->length - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}

	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_z + 1][pos_y]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x][pos_z]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][route_x->length - 1]);
}

inline void LocalSearch::ConcatRouteXRelocate2IntraBack()
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_y; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	//pathX.push_back(cl_x);
	//pathX.push_back(cl_z);
	for (int i = pos_x; i <= pos_z; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}

	for (int i = pos_y + 1; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_z + 1; i <= route_x->length - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_y], route_x->subseq[pos_x][pos_z]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][pos_x - 1]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_z + 1][route_x->length - 1]);
}

inline void LocalSearch::EvalRelocate2()
{
	double temp = -data->time_cost[prev_x][cl_x] - data->time_cost[cl_z][next_z] + data->time_cost[prev_x][next_z] -
		data->time_cost[cl_y][cl_w];
	curr_delta1 = temp + data->time_cost[cl_y][cl_x] + data->time_cost[cl_z][cl_w];
	curr_delta2 = temp + data->time_cost[cl_y][cl_z] + data->time_cost[cl_x][cl_w];
	curr_delta=0;

	if (route_x->is_feasible && route_y->is_feasible && curr_delta1 >= EPSILON && curr_delta2 >= EPSILON)
		return;

#ifdef TW
	ConcatRouteXRelocate2Inter();
	ConcatRouteYRelocate2Inter();

#ifdef ROBUST_TW
	curr_delta1 += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolation);
	curr_delta2 += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolationRev);
#else
	temp = data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
	curr_delta1 += temp + data->penalty_tw * (-route_y->tw_violation + subseq_y.timewarp);
	curr_delta2 += temp + data->penalty_tw * (-route_y->tw_violation + subseq_y_rev.timewarp);
#endif
#endif

#ifdef ROBUST_DEMAND
	ComputeRobustDemandRouteXRelocate2();
	ComputeRobustDemandRouteYRelocate2();
	temp = data->penalty_load * (routeXRobustCapViolation + routeYRobustCapViolation -
			route_x->robust_cap_violation - route_y->robust_cap_violation);
	curr_delta1 += temp;
	curr_delta2 += temp;
#else
	double tempLoadDelta = -data->client[cl_x].demand - data->client[cl_z].demand;
	temp = data->penalty_load * (-route_x->cap_violation + ComputeCapViolation(route_x->load.back() + tempLoadDelta)) +
		data->penalty_load * (-route_y->cap_violation + ComputeCapViolation(route_y->load.back() - tempLoadDelta));
	curr_delta1 += temp;
	curr_delta2 += temp;
#endif

	if (curr_delta1 < curr_delta2)
	{
		curr_delta = curr_delta1;
		move_type = 0;
	}
	else
	{
		curr_delta = curr_delta2;
		move_type = 1;
	}

#ifdef USEMMD
	if (curr_delta + EPSILON < inter_smd->cost)
	{
		inter_smd->status = IMPROVED;
		inter_smd->cost = curr_delta;
		inter_smd->r1_pos1 = pos_x;
		inter_smd->r2_pos1 = pos_y;
		inter_smd->move_type = move_type;
	}
#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
		bestMoveType = move_type;
	}
#endif
}

inline void LocalSearch::EvalOrOpt()
{
	curr_delta = -data->time_cost[prev_x][cl_x] - data->time_cost[cl_z][next_z] + data->time_cost[prev_x][next_z] -
		data->time_cost[cl_y][cl_w] + data->time_cost[cl_y][cl_x] + data->time_cost[cl_z][cl_w];

#ifdef TW
	if (route_x->is_time_feasible && curr_delta >= EPSILON)
		return;

	if (pos_x < pos_y)
		ConcatRouteXRelocate2IntraFront();
	else
		ConcatRouteXRelocate2IntraBack();

#ifdef ROBUST_TW
	curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation);
#else
	curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
#endif
#endif

#ifdef USEMMD
	if (curr_delta + EPSILON < intra_smd->cost)
	{
		intra_smd->status = IMPROVED;
		intra_smd->cost = curr_delta;
		intra_smd->r1_pos1 = pos_x;
		intra_smd->r2_pos1 = pos_y;
	}
#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
	}
#endif
}

inline void LocalSearch::ExecRelocate2()
{
	if (best_route_x != best_route_y)
	{
		for (int i = best_pos_x; i < best_route_x->length - 2; i++)
		{
			std::swap(curr_solution->route[best_route_x->cour][i + 1], curr_solution->route[best_route_x->cour][i + 2]);
			std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);
		}
		int cust1, cust2;
		if (bestMoveType == 0)
		{
			cust2 = curr_solution->route[best_route_x->cour].back();
			curr_solution->route[best_route_x->cour].pop_back();
			cust1 = curr_solution->route[best_route_x->cour].back();
			curr_solution->route[best_route_x->cour].pop_back();
		}
		else
		{
			cust1 = curr_solution->route[best_route_x->cour].back();
			curr_solution->route[best_route_x->cour].pop_back();
			cust2 = curr_solution->route[best_route_x->cour].back();
			curr_solution->route[best_route_x->cour].pop_back();
		}

		best_route_x->Update();

		curr_solution->route[best_route_y->cour].push_back(cust1);
		curr_solution->route[best_route_y->cour].push_back(cust2);
		for (int i = best_route_y->length - 1; i > best_pos_y; i--)
		{
			std::swap(curr_solution->route[best_route_y->cour][i], curr_solution->route[best_route_y->cour][i + 1]);
			std::swap(curr_solution->route[best_route_y->cour][i + 1], curr_solution->route[best_route_y->cour][i + 2]);
		}
		best_route_y->Update();
	}

	else
	{
		if (best_pos_x < best_pos_y)
			for (int i = best_pos_x; i <= best_pos_y - 2; i++)
			{
				std::swap(curr_solution->route[best_route_x->cour][i + 1], curr_solution->route[best_route_x->cour][i + 2]);
				std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);
			}
		else
			for (int i = best_pos_x; i >= best_pos_y + 2; i--)
			{
				std::swap(curr_solution->route[best_route_x->cour][i - 1], curr_solution->route[best_route_x->cour][i]);
				std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);
			}
		best_route_x->Update();
	}

	UpdateClientRank();

	// std::cout << *curr_solution << "\n";
}

void LocalSearch::SearchRelocate2()
{
	//     improved = false;
	//     best_delta = 0;
	//     for (cl_x = 1; cl_x <= data->nb_clients; cl_x++)
	//     {
	//         pos_x = clientRank[cl_x].second;
	//         routeX = clientRank[cl_x].first;
	//         pos_z = pos_x + 1;
	//         cl_z = curr_solution->route[rx_idx][pos_z];
	//         if (cl_z == 0)
	//             continue;
	// 
	//         for (int cy = 0; cy < data->correlated_vertices[cl_x].size(); cy++)
	//         {
	//             cl_y = data->correlated_vertices[cl_x][cy];
	//             pos_y = clientRank[cl_y].second;
	//             pos_w = pos_y + 1;
	//             routeY = clientRank[cl_y].first;
	//             cl_w = curr_solution->route[ry_idx][pos_w];
	//             prev_x = curr_solution->route[rx_idx][pos_x - 1];
	//             prev_y = curr_solution->route[ry_idx][pos_y - (pos_y > 0)];
	//             next_z = curr_solution->route[rx_idx][pos_z + 1];
	//             // next_w = curr_solution->route[ry_idx][pos_w + 1];
	// 
	//             // if(routeX != routeY) continue;
	// 
	//             // SetLocalVariables();
	//             EvalRelocate2();
	//         }
	//     }
	// 
	//     // for (int r1 = 0; r1 < ads->routeADS.size(); r1++)
	//     // {
	//     //     for (int r2 = 0; r2 < ads->routeADS.size(); r2++)
	//     //     {
	//     //         if (routeData[r1].path.size() == 2 || r1 == r2)
	//     //             continue;
	// 
	//     //         for (int i = 1; i < routeData[r1].path.size() - 1; i++)
	//     //         {
	//     //             for (int j = 0; j < routeData[r2].path.size() - 1; j++)
	//     //             {
	//     //                 cl_x = routeData[r1].path[i];
	//     //                 cl_y = routeData[r2].path[j];
	//     //                 pos_z = pos_x = clientRank[cl_x].second;
	//     //                 pos_w = pos_y = clientRank[cl_y].second;
	//     //                 pos_z++, pos_w++;
	//     //                 routeX = clientRank[cl_x].first;
	//     //                 routeY = clientRank[cl_y].first;
	//     //                 cl_z = curr_solution->route[rx_idx][pos_z];
	//     //                 cl_w = curr_solution->route[ry_idx][pos_w];
	//     //                 prev_x = curr_solution->route[rx_idx][pos_x - 1];
	//     //                 prev_y = curr_solution->route[ry_idx][pos_y - (pos_y > 0)];
	//     //                 EvalRelocate1();
	//     //             }
	//     //         }
	//     //     }
	//     // }
	// 
	//     if (improved)
	//     {
	//         // std::cout << best_delta << "\n";
	//         // std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";
	//         // std::cout << *curr_solution << "\n";
	//         ExecRelocate2();
	//     }
}

void LocalSearch::SearchRelocate2Complete()
{
	stats->inter_neighb_explore[4].SetStart();
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		route_x = ads->routeADS[rx_idx];
		if (route_x->length <= 3)
			continue;

		for (ry_idx = 0; ry_idx < ads->routeADS.size(); ry_idx++)
		{
			route_y = ads->routeADS[ry_idx];
			rx_ry_feasible = route_x->is_feasible && route_y->is_feasible;
			rx_ry_improved = false;

			//if (rx_idx == ry_idx) continue;
			if (rx_idx == ry_idx 
				|| (rx_ry_feasible && (route_y->load.back() + route_x->min_pair_demand > data->capacity + EPSILON))
				|| (!curr_solution->GetInterNeighbStatus(rx_idx, ry_idx, 4)))
				continue;


			stats->inter_mmd_search[4].SetStart();
#ifdef USEMMD
			inter_smd = ads->getSMD(rx_idx, ry_idx, 4);	
#endif
			stats->inter_mmd_search[4].SetEnd();

			start_search = 0;
			end_search  = route_y->length - 2;

			if (!use_mmd || inter_smd->status == UNVISITED)
			{
#ifdef USEMMD
				inter_smd->status = UNIMPROVED;
#endif

				for (pos_x = 1; pos_x < route_x->length - 2; pos_x++)
				{
					pos_z = pos_x + 1;
					cl_x = curr_solution->route[rx_idx][pos_x];
					cl_z = curr_solution->route[rx_idx][pos_z];
					if (route_y->length <= 2 || rx_ry_feasible && (route_y->load.back() + data->client[cl_x].demand + data->client[cl_z].demand > data->capacity))
						continue;
					prev_x = curr_solution->route[rx_idx][pos_x - 1];
					next_z = curr_solution->route[rx_idx][pos_z + 1];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_x];
						end_search = route_y->end_search[cl_x];
						// std::cout << "\n";
						// std::cout << start_search << " " << 0 << "\n";
						// std::cout << end_search << " " << route_y->length - 2 << "\n";
						// std::cout << "\n";
					}
#endif

					// for (pos_y = 0; pos_y < route_y->length - 1; pos_y++)
					for (pos_y = start_search; pos_y <= end_search; pos_y++)
					{
						pos_w = pos_y + 1;
						cl_y = curr_solution->route[ry_idx][pos_y];
						cl_w = curr_solution->route[ry_idx][pos_w];
						// prev_y = curr_solution->route[ry_idx][pos_y - (pos_y > 0)];

						EvalRelocate2();
						rx_ry_improved |= (curr_delta + EPSILON < 0);
					}
				}
			}
			else
			{
				stats->inter_mmd_effective[4]++;
			}
#ifdef USEMMD		
			if (inter_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = inter_smd->cost;
				best_pos_x = inter_smd->r1_pos1;
				best_pos_y = inter_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
				bestMoveType = inter_smd->move_type;
			}
			rx_ry_improved |= (inter_smd->cost + EPSILON < 0);
#endif

			if (!rx_ry_improved)
				curr_solution->SetInterNeighbStatus(rx_idx, ry_idx, 4, false);
		}
	}

	if (improved)
	{
		stats->neighb_improv[1]++;
		// std::cout << *curr_solution << "\n";
		// std::cout << best_delta << "\n";
		// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x->cour << " " << best_route_y->cour << "\n";

		curr_solution->SetInterNeighbStatusAllTrue(best_route_x->cour, best_route_y->cour);
		ExecRelocate2();
	}
	stats->inter_neighb_explore[4].SetEnd();
}

/** ********************* SWAP (2,1) ***************************** */

inline void LocalSearch::ComputeRobustDemandRouteXSwap21()
{
	routeXRobustDemand = route_x->load.back() - data->client[cl_x].demand -
		data->client[cl_z].demand + data->client[cl_y].demand;
	std::vector<double> sorted_deviations(route_x->length - 3);

	int j = 0;
	for (int i = 1; i <= pos_x - 1; i++, j++)
	{
		sorted_deviations[j] = -data->client[curr_solution->route[rx_idx][i]].demand_deviation;
	}
	for (int i = pos_x + 2; i <= route_x->length - 2; i++, j++)
	{
		sorted_deviations[j] = -data->client[curr_solution->route[rx_idx][i]].demand_deviation;
	}
	sorted_deviations[j] = -data->client[cl_y].demand_deviation;

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeXRobustDemand -= sorted_deviations[i];
		}
	}
	routeXRobustCapViolation = std::max(0.0, routeXRobustDemand - data->capacity);
}

inline void LocalSearch::ComputeRobustDemandRouteYSwap21()
{
	routeYRobustDemand = route_y->load.back() + data->client[cl_x].demand +
		data->client[cl_z].demand - data->client[cl_y].demand;
	std::vector<double> sorted_deviations(route_y->length - 1);

	int j = 0;
	for (int i = 1; i <= pos_y - 1; i++, j++)
	{
		sorted_deviations[j] = -data->client[curr_solution->route[ry_idx][i]].demand_deviation;
	}
	for (int i = pos_y + 1; i <= route_y->length - 2; i++, j++)
	{
		sorted_deviations[j] = -data->client[curr_solution->route[ry_idx][i]].demand_deviation;
	}
	sorted_deviations[j++] = -data->client[cl_x].demand_deviation;
	sorted_deviations[j++] = -data->client[cl_z].demand_deviation;

	if (sorted_deviations.size() > data->budget_demand)
	{
		std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin() + data->budget_demand, sorted_deviations.end());

		for (int i = 0; i < data->budget_demand; i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}
	else
	{
		for (int i = 0; i < sorted_deviations.size(); i++)
		{
			routeYRobustDemand -= sorted_deviations[i];
		}
	}

	routeYRobustCapViolation = std::max(0.0, routeYRobustDemand - data->capacity);
}

inline void LocalSearch::ConcatRouteXSwap21Inter()
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	for (int i = pos_x + 2; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_y->subseq[pos_y][pos_y]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 2][route_x->length - 1]);
}

inline void LocalSearch::ConcatRouteYSwap21Inter()
{
#ifdef ROBUST_TW
	routeYRobustTWViolation = 0;
	pathY.clear();
	for (int i = 0; i <= pos_y - 1; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	pathY.push_back(cl_x);
	pathY.push_back(cl_z);
	for (int i = pos_y + 1; i < route_y->length; i++)
	{
		pathY.push_back(curr_solution->route[ry_idx][i]);
	}
	ComputePathEarliestStart(pathY, routeYRobustTWViolation);
	routeYRobustTWViolationRev = 0;
	std::swap(pathY[pos_y+1], pathY[pos_y]);
	ComputePathEarliestStart(pathY, routeYRobustTWViolationRev);
#endif

	subseq_y = Subsequence::Concatenate(data, route_y->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x + 1]);
	subseq_y = Subsequence::Concatenate(data, subseq_y, route_y->subseq[pos_y + 1][route_y->length - 1]);
	subseq_y_rev = Subsequence::Concatenate(data, route_y->subseq[0][pos_y - 1], route_x->subseq[pos_x + 1][pos_x]);
	subseq_y_rev = Subsequence::Concatenate(data, subseq_y_rev, route_y->subseq[pos_y + 1][route_y->length - 1]);
}

inline void LocalSearch::ConcatRouteXSwap21IntraFront() // when pos_x < pos_y
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	for (int i = pos_x + 2; i <= pos_y - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_x);
	pathX.push_back(cl_z);
	for (int i = pos_y + 1; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	if (pos_z + 1 == pos_y)
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_y][pos_y]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x][pos_x + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][route_x->length - 1]);
	}
	else
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x - 1], route_x->subseq[pos_y][pos_y]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 2][pos_y - 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x][pos_x + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][route_x->length - 1]);
	}
}
inline void LocalSearch::ConcatRouteXSwap21IntraBack() // when pos_x > pos_y
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_y - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_x);
	pathX.push_back(cl_z);
	for (int i = pos_y + 1; i <= pos_x - 1; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	pathX.push_back(cl_y);
	for (int i = pos_x + 2; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	if (pos_y + 1 == pos_x)
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y][pos_y]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 2][route_x->length - 1]);
	}
	else
	{
		subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_y - 1], route_x->subseq[pos_x][pos_x + 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][pos_x - 1]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y][pos_y]);
		subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_x + 2][route_x->length - 1]);
	}
}

inline void LocalSearch::EvalSwap21()
{
//	if (rx_idx != ry_idx)
//	{
	double temp = data->time_cost[prev_x][cl_y] + data->time_cost[cl_y][next_z] 
							- data->time_cost[prev_y][cl_y] - data->time_cost[cl_y][cl_w]
							+ temp_delta;
	curr_delta1 = temp + data->time_cost[prev_y][cl_x] + data->time_cost[cl_z][cl_w];
	curr_delta2 = temp + data->time_cost[prev_y][cl_z] + data->time_cost[cl_x][cl_w];
	curr_delta = 0;

	if (route_x->is_feasible && route_y->is_feasible && curr_delta1 >= EPSILON && curr_delta2 >= EPSILON)
		return;

#ifdef TW
	ConcatRouteXSwap21Inter();
	ConcatRouteYSwap21Inter();

#ifdef ROBUST_TW
	curr_delta1 += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolation);
	curr_delta2 += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation) +
		data->penalty_tw * (-route_y->robust_tw_violation + routeYRobustTWViolationRev);
#else
	temp = data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
	curr_delta1 += temp + data->penalty_tw * (-route_y->tw_violation + subseq_y.timewarp);
	curr_delta2 += temp + data->penalty_tw * (-route_y->tw_violation + subseq_y_rev.timewarp);
#endif
#endif

#ifdef ROBUST_DEMAND
	ComputeRobustDemandRouteXSwap21();
	ComputeRobustDemandRouteYSwap21();
	temp = data->penalty_load * (routeXRobustCapViolation + routeYRobustCapViolation -
			route_x->robust_cap_violation - route_y->robust_cap_violation);
	curr_delta1 += temp;
	curr_delta2 += temp;
#else
	double tempLoadDelta = -data->client[cl_x].demand - data->client[cl_z].demand + data->client[cl_y].demand;
	temp = data->penalty_load * (-route_x->cap_violation + ComputeCapViolation(route_x->load.back() + tempLoadDelta)) +
		data->penalty_load * (-route_y->cap_violation + ComputeCapViolation(route_y->load.back() - tempLoadDelta));
	curr_delta1 += temp;
	curr_delta2 += temp;
#endif

	curr_delta = curr_delta1;
	move_type = 0;

	if (curr_delta2 + EPSILON < curr_delta)
	{
		curr_delta = curr_delta2;
		move_type = 1;
	}

	//if (curr_delta1 + EPSILON < curr_delta2)
	//{
	//	curr_delta = curr_delta1;
	//	move_type = 0;
	//}
	//else
	//{
	//	curr_delta = curr_delta2;
	//	move_type = 1;
	//}

	//if (std::abs(curr_delta1 - curr_delta2) <= EPSILON)
	//{
	//	curr_delta = curr_delta1;
	//	move_type = 0;
	//}

#ifdef USEMMD
	if (curr_delta + EPSILON < inter_smd->cost)
	{
		inter_smd->status = IMPROVED;
		inter_smd->cost = curr_delta;
		inter_smd->r1_pos1 = pos_x;
		inter_smd->r2_pos1 = pos_y;
		inter_smd->move_type = move_type;
	}
#else

	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
		bestMoveType = move_type;
	}
#endif

//	}
//	else if (cl_y != cl_z)
//	{
//		if (pos_z + 1 == pos_y)
//		{
//			curr_delta = -data->time_cost[prev_x][cl_x] - data->time_cost[cl_z][cl_y] + data->time_cost[prev_x][cl_y] + data->time_cost[cl_y][cl_x] -
//				data->time_cost[cl_y][cl_w] + data->time_cost[cl_z][cl_w];
//		}
//		else if (pos_y + 1 == pos_x)
//		{
//			curr_delta = -data->time_cost[prev_y][cl_y] - data->time_cost[cl_y][cl_x] + data->time_cost[cl_z][cl_y] + data->time_cost[cl_y][next_z] -
//				data->time_cost[cl_z][next_z] + data->time_cost[prev_y][cl_x];
//		}
//		else
//		{
//			curr_delta = -data->time_cost[prev_x][cl_x] - data->time_cost[cl_z][next_z] + data->time_cost[prev_x][cl_y] + data->time_cost[cl_y][next_z] +
//				-data->time_cost[prev_y][cl_y] - data->time_cost[cl_y][cl_w] + data->time_cost[prev_y][cl_x] + data->time_cost[cl_z][cl_w];
//		}
//
//#ifdef TW
//		if (route_x->is_time_feasible && curr_delta >= EPSILON)
//			return;
//
//		if (pos_x < pos_y)
//		{
//			ConcatRouteXSwap21IntraFront();
//		}
//		else
//		{
//			ConcatRouteXSwap21IntraBack();
//		}
//
//#ifdef ROBUST_TW
//		curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation);
//#else
//		curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
//#endif
//#endif
//
//		if (curr_delta + EPSILON < best_delta)
//		{
//			improved = true;
//			best_delta = curr_delta;
//			best_pos_x = pos_x;
//			best_pos_y = pos_y;
//			best_route_x = route_x;
//			best_route_y = route_y;
//			// std::cout << "\n\n\n\n";
//			// std::cout << *curr_solution << "\n";
//			// std::cout << "BestDelta: " << best_delta << "\n";
//			// std::cout << "BestPosX: " << best_pos_x << '\n';
//			// std::cout << "BEstPoSY: " << best_pos_y << '\n';
//			// std::cout << "RouteX: " << best_route_x << '\n';
//			// std::cout << "RouteY: " << best_route_y << '\n';
//			// std::cout << "\n\n\n\n";
//		}
//	}
}

inline void LocalSearch::ExecSwap21()
{
	if (best_route_x != best_route_y)
	{
		int cust1, cust2;
		if (bestMoveType == 0)
		{
			cust1 = curr_solution->route[best_route_x->cour][best_pos_x];
			cust2 = curr_solution->route[best_route_x->cour][best_pos_x + 1];
		}
		else
		{
			cust2 = curr_solution->route[best_route_x->cour][best_pos_x];
			cust1 = curr_solution->route[best_route_x->cour][best_pos_x + 1];
		}
		curr_solution->route[best_route_x->cour][best_pos_x] = curr_solution->route[best_route_y->cour][best_pos_y];

		for (int i = best_pos_x + 1; i < best_route_x->length - 1; i++)
			std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);

		curr_solution->route[best_route_x->cour].pop_back();

		curr_solution->route[best_route_y->cour][best_pos_y] = cust1;
		curr_solution->route[best_route_y->cour].push_back(cust2);

		for (int i = best_route_y->length - 1; i > best_pos_y; i--)
			std::swap(curr_solution->route[best_route_y->cour][i], curr_solution->route[best_route_y->cour][i + 1]);

		best_route_x->Update();
		best_route_y->Update();
	}

	else
	{
		if (best_pos_x < best_pos_y)
		{
			std::swap(curr_solution->route[best_route_x->cour][best_pos_x], curr_solution->route[best_route_x->cour][best_pos_y]);
			for (int i = best_pos_x + 1; i <= best_pos_y - 1; i++)
			{
				std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);
			}
		}
		else
		{
			std::swap(curr_solution->route[best_route_x->cour][best_pos_x + 1], curr_solution->route[best_route_x->cour][best_pos_y]);
			for (int i = best_pos_x; i >= best_pos_y + 1; i--)
			{
				std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i - 1]);
			}
		}

		best_route_x->Update();
	}
	UpdateClientRank();
}

inline void LocalSearch::SearchSwap21()
{
	//     improved = false;
	//     best_delta = 0;
	//     for (cl_x = 1; cl_x <= data->nb_clients; cl_x++)
	//     {
	//         for (int cy = 0; cy < data->correlated_vertices[cl_x].size(); cy++)
	//         {
	//             prev_y = data->correlated_vertices[cl_x][cy];
	//             routeY = clientRank[prev_y].first;
	//             pos_y = clientRank[prev_y].second + 1;
	//             cl_y = curr_solution->route[ry_idx][pos_y];
	// 
	//             if (cl_y == cl_x)
	//                 continue;
	// 
	//             if (pos_y >= route_y->length - 3 || route_y->length < 4)
	//                 continue;
	// 
	//             pos_x = clientRank[cl_x].second;
	//             pos_z = pos_x + 1;
	//             pos_w = pos_y + 1;
	//             routeX = clientRank[cl_x].first;
	//             if (pos_x >= route_x->length - 3 || route_x->length < 4)
	//                 continue;
	// 
	//             cl_z = curr_solution->route[rx_idx][pos_z];
	//             cl_w = curr_solution->route[ry_idx][pos_w];
	//             prev_x = curr_solution->route[rx_idx][pos_x - 1];
	//             next_w = curr_solution->route[ry_idx][pos_w + 1];
	//             next_z = curr_solution->route[rx_idx][pos_z + 1];
	// 
	//             // if (routeX != routeY)
	//             //     continue;
	// 
	//             // if (pos_x < pos_y)
	//             // {
	//             //     if(pos_x + 2 == pos_y)
	//             //     {
	//             //         std::cout << "A\n";
	//             //     }
	//             // }
	// 
	//             // SetLocalVariables();
	//             EvalSwap21();
	//         }
	//     }
	//     if (improved)
	//     {
	//         // std::cout << best_delta << "\n";
	//         // std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";
	//         // std::cout << *curr_solution;
	//         ExecSwap21();
	//     }
}

void LocalSearch::SearchSwap21Complete()
{
	stats->inter_neighb_explore[5].SetStart();
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		route_x = ads->routeADS[rx_idx];
		if (route_x->length <= 3)
			continue;

		for (ry_idx = 0; ry_idx < ads->routeADS.size(); ry_idx++)
		{
			route_y = ads->routeADS[ry_idx];
			rx_ry_feasible = route_x->is_feasible && route_y->is_feasible;
			rx_ry_improved  = false;
			if (rx_idx == ry_idx || route_y->length <= 2 || route_x->length <= 3
						|| (rx_ry_feasible && (route_x->load.back() - route_x->max_pair_demand + route_y->min_demand > data->capacity + EPSILON))
				 || (!curr_solution->GetInterNeighbStatus(rx_idx, ry_idx, 5)))
				continue;

			stats->inter_mmd_search[5].SetStart();
#ifdef USEMMD
			inter_smd = ads->getSMD(rx_idx, ry_idx, 5);	
#endif
			stats->inter_mmd_search[5].SetEnd();

			start_search = 0;
			end_search  = route_y->length - 2;

			if (!use_mmd || inter_smd->status == UNVISITED)
			{
#ifdef USEMMD
				inter_smd->status = UNIMPROVED;
#endif

				for (pos_x = 1; pos_x < route_x->length - 2; pos_x++)
				{
					pos_z = pos_x + 1;
					cl_x = curr_solution->route[rx_idx][pos_x];
					cl_z = curr_solution->route[rx_idx][pos_z];
					if (rx_ry_feasible && (route_y->load.back() - route_y->max_demand
								+ data->client[cl_x].demand 
								+ data->client[cl_z].demand > data->capacity))
						continue;

					prev_x = curr_solution->route[rx_idx][pos_x - 1];
					next_z = curr_solution->route[rx_idx][pos_z + 1];

					temp_delta = - data->time_cost[prev_x][cl_x] - data->time_cost[cl_z][next_z];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_x];
						end_search = route_y->end_search[cl_x];
						// std::cout << "\n";
						// std::cout << start_search << " " << 0 << "\n";
						// std::cout << end_search << " " << route_y->length - 2 << "\n";
						// std::cout << "\n";
					}
#endif

					for (pos_y = start_search + 1; pos_y <= end_search; pos_y++)
					{
						pos_w = pos_y + 1;
						cl_y = curr_solution->route[ry_idx][pos_y];
						cl_w = curr_solution->route[ry_idx][pos_w];
						prev_y = curr_solution->route[ry_idx][pos_y - 1];
						// next_w = curr_solution->route[ry_idx][pos_w + 1];
						EvalSwap21();
						rx_ry_improved |= (curr_delta + EPSILON < 0);
					}
				}
			}
			else
			{
				stats->inter_mmd_effective[5]++;
			}
#ifdef USEMMD		
			if (inter_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = inter_smd->cost;
				best_pos_x = inter_smd->r1_pos1;
				best_pos_y = inter_smd->r2_pos1;
				bestMoveType = inter_smd->move_type;
				best_route_x = route_x;
				best_route_y = route_y;
			}
			rx_ry_improved |= (inter_smd->cost + EPSILON < 0);
#endif

			if (!rx_ry_improved)
				curr_solution->SetInterNeighbStatus(rx_idx, ry_idx, 5, false);
		}
	}

	if (improved)
	{
		// std::cout << best_delta << " " << best_pos_x << " " << best_pos_y << " " << bestMoveType << "\n";
		stats->neighb_improv[3]++;
		// std::cout << *curr_solution << "\n";
		// std::cout << best_delta << "\n";
		// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";

		curr_solution->SetInterNeighbStatusAllTrue(best_route_x->cour, best_route_y->cour);
		ExecSwap21();
	}
	stats->inter_neighb_explore[5].SetEnd();
}

inline void LocalSearch::ConcatRouteXTwoOpt()
{
#ifdef ROBUST_TW
	routeXRobustTWViolation = 0;
	pathX.clear();
	for (int i = 0; i <= pos_x; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_y; i > pos_x; i--)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	for (int i = pos_y + 1; i < route_x->length; i++)
	{
		pathX.push_back(curr_solution->route[rx_idx][i]);
	}
	ComputePathEarliestStart(pathX, routeXRobustTWViolation);
#endif

	subseq_x = Subsequence::Concatenate(data, route_x->subseq[0][pos_x], route_x->subseq[pos_y][pos_x + 1]);
	subseq_x = Subsequence::Concatenate(data, subseq_x, route_x->subseq[pos_y + 1][route_x->length - 1]);
}

inline void LocalSearch::EvalTwoOpt()
{
	curr_delta = -data->time_cost[cl_x][cl_z] - data->time_cost[cl_y][cl_w] +
		data->time_cost[cl_x][cl_y] + data->time_cost[cl_z][cl_w];

#ifdef TW
	if (route_x->is_time_feasible && curr_delta >= EPSILON)
		return;

	ConcatRouteXTwoOpt();

#ifdef ROBUST_TW
	curr_delta += data->penalty_tw * (-route_x->robust_tw_violation + routeXRobustTWViolation);
#else
	curr_delta += data->penalty_tw * (-route_x->tw_violation + subseq_x.timewarp);
#endif
#endif

#ifdef USEMMD
	if (curr_delta + EPSILON < intra_smd->cost)
	{
		intra_smd->status = IMPROVED;
		intra_smd->cost = curr_delta;
		intra_smd->r1_pos1 = pos_x;
		intra_smd->r2_pos1 = pos_y;
	}
#else
	if (curr_delta + EPSILON < best_delta)
	{
		improved = true;
		best_delta = curr_delta;
		best_pos_x = pos_x;
		best_pos_y = pos_y;
		best_route_x = route_x;
		best_route_y = route_y;
	}
#endif
}
inline void LocalSearch::ExecTwoOpt()
{
	std::reverse(curr_solution->route[best_route_x->cour].begin() + best_pos_x + 1, curr_solution->route[best_route_x->cour].begin() + best_pos_y + 1);
	best_route_x->Update();
	UpdateClientRank();
}

inline void LocalSearch::SearchTwoOpt()
{
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		bool changed = false;
		ry_idx = rx_idx;
//		improved = false;
		best_delta = 0;
		route_x = ads->routeADS[rx_idx];
		route_y = ads->routeADS[ry_idx];

		int counter = 0;

		while (counter < 100)
		{
			if (route_x->length <= 6
					|| !curr_solution->GetIntraNeighbStatus(rx_idx, 4))
				break;

			improved = false;

#ifdef USEMMD
			intra_smd = ads->getSMD(rx_idx, 4);
#endif

			if (!use_mmd || intra_smd->status == UNVISITED)
			{
#ifdef USEMMD
				intra_smd->status = UNIMPROVED;
#endif
				for (pos_x = 1; pos_x < route_x->length - 2; pos_x++)
				{
					pos_z = pos_x + 1;
					cl_x = curr_solution->route[rx_idx][pos_x];
					cl_z = curr_solution->route[rx_idx][pos_z];
					prev_x = curr_solution->route[rx_idx][pos_x - 1];

					for (pos_y = pos_x + 1; pos_y < route_x->length - 1; pos_y++)
					{
						pos_w = pos_y + 1;
						cl_y = curr_solution->route[rx_idx][pos_y];
						cl_w = curr_solution->route[rx_idx][pos_w];
						prev_y = curr_solution->route[rx_idx][pos_y - (pos_y > 0)];
						EvalTwoOpt();
					}
				}
			}
#ifdef USEMMD		
			if (intra_smd->cost + EPSILON < best_delta)
			{
				// std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << route_x->random_cost << "\n";
				// std::cout << "found " << std::setprecision(std::numeric_limits<double>::max_digits10) << route_x->random_cost << " " << route_x->length << "\n";
				improved = true;
				best_delta = intra_smd->cost;
				best_pos_x = intra_smd->r1_pos1;
				best_pos_y = intra_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
			}
#endif

			if (improved)
			{
				stats->neighb_improv[10]++;
				// std::cout << best_delta << "\n";
				// std::cout << *curr_solution << "\n";
				// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";
				// double cost_check1 = best_route_x->distance;
				// cost_check1 += data->penalty_tw * best_route_x->robust_tw_violation;
				// cost_check1 += data->penalty_load * best_route_x->robust_cap_violation;

				ExecTwoOpt();
				//curr_solution->SetInterNeighbStatusAllTrue(rx_idx);
				//curr_solution->SetIntraNeighbStatusAllTrue(rx_idx);

				// double cost_check2 = best_route_x->distance;
				// cost_check2 += data->penalty_tw * best_route_x->robust_tw_violation;
				// cost_check2 += data->penalty_load * best_route_x->robust_cap_violation;

				// std::cout << cost_check1 - cost_check2 << " " << best_delta << "\n";
			}

			else
			{
				counter = 100;
				break;
			}

			counter++;
		}
	}
}

void LocalSearch::SearchExchange()
{
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		ry_idx = rx_idx;
		//improved = false;
		best_delta = 0;
		route_x = ads->routeADS[rx_idx];
		route_y = ads->routeADS[ry_idx];

		int counter = 0;

		while (counter < 100)
		{
			if (route_x->length <= 3
					|| !curr_solution->GetIntraNeighbStatus(rx_idx, 0))
				break;

			improved = false;

#ifdef USEMMD
			intra_smd = ads->getSMD(rx_idx, 0);	
#endif

			if (!use_mmd || intra_smd->status == UNVISITED)
			{
#ifdef USEMMD
				intra_smd->status = UNIMPROVED;
#endif
				for (pos_x = 1; pos_x < route_x->length - 2; pos_x++)
				{
					pos_z = pos_x + 1;
					cl_x = curr_solution->route[rx_idx][pos_x];
					cl_z = curr_solution->route[rx_idx][pos_z];
					prev_x = curr_solution->route[rx_idx][pos_x - 1];

					for (pos_y = pos_x + 1; pos_y < route_x->length - 1; pos_y++)
					{
						pos_w = pos_y + 1;
						cl_y = curr_solution->route[rx_idx][pos_y];
						cl_w = curr_solution->route[rx_idx][pos_w];
						prev_y = curr_solution->route[rx_idx][pos_y - (pos_y > 0)];
						EvalExchange();
					}
				}
			}
#ifdef USEMMD		
			if (intra_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = intra_smd->cost;
				best_pos_x = intra_smd->r1_pos1;
				best_pos_y = intra_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
			}
#endif

			if (improved)
			{
				stats->neighb_improv[6]++;
				// std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << best_route_x->random_cost << " " << best_route_x->length << "\n";
				// std::cout << best_delta << "\n";
				// std::cout << *curr_solution << "\n";
				// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";

				//curr_solution->SetInterNeighbStatusAllTrue(rx_idx);
				//curr_solution->SetIntraNeighbStatusAllTrue(rx_idx);
				ExecSwap11();
			}

			else
			{
				counter = 100;
				break;
			}

			counter++;
		}
	}

}

void LocalSearch::SearchReinsertion()
{
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
	//	improved = false;
		best_delta = 0;
		ry_idx = rx_idx;
		route_x = ads->routeADS[rx_idx];
		route_y = ads->routeADS[ry_idx];

		int counter = 0;

		while (counter < 100)
		{
			if (route_x->length <= 4
					|| !curr_solution->GetIntraNeighbStatus(rx_idx, 1))
				break;

			improved = false;

#ifdef USEMMD
			intra_smd = ads->getSMD(rx_idx, 1);	
#endif

			start_search = 0;
			end_search  = route_x->length - 2;
			rx_ry_feasible = route_x->is_feasible;

			if (!use_mmd || intra_smd->status == UNVISITED)
			{
#ifdef USEMMD
				intra_smd->status = UNIMPROVED;
#endif
				for (pos_x = 1; pos_x < route_x->length - 1; pos_x++)
				{
					cl_x = curr_solution->route[rx_idx][pos_x];
					pos_z = pos_x + 1;
					cl_z = curr_solution->route[rx_idx][pos_z];
					prev_x = curr_solution->route[rx_idx][pos_x - 1];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_x];
						end_search = route_y->end_search[cl_x];
						// std::cout << "\n";
						// std::cout << start_search << " " << 0 << "\n";
						// std::cout << end_search << " " << route_y->length - 2 << "\n";
						// std::cout << "\n";
					}
#endif

					for (pos_y = start_search; pos_y <= end_search; pos_y++)
					{
						pos_w = pos_y + 1;
						if (pos_x == pos_w || pos_y == pos_x)
							continue;
						cl_y = curr_solution->route[rx_idx][pos_y];
						cl_w = curr_solution->route[rx_idx][pos_w];
						// prev_y = curr_solution->route[rx_idx][pos_y - 1];
						EvalReinsertion();
					}
				}
			}

#ifdef USEMMD		
			if (intra_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = intra_smd->cost;
				best_pos_x = intra_smd->r1_pos1;
				best_pos_y = intra_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
			}
#endif

			if (improved)
			{
				stats->neighb_improv[7]++;
				ExecRelocate1();
				//curr_solution->SetInterNeighbStatusAllTrue(rx_idx);
				//curr_solution->SetIntraNeighbStatusAllTrue(rx_idx);
				// 	double cost_check2 = best_route_x->distance;
				// 	cost_check2 += data->penalty_tw * best_route_x->robust_tw_violation;
				// 	cost_check2 += data->penalty_load * best_route_x->robust_cap_violation;

				//   std::cout << *curr_solution << "\n";
				// 	std::cout << cost_check1 - cost_check2 << " " << best_delta << "\n";
			}

			else
			{
				counter = 100;
				break;
			}
			
			counter++;
		}
	}
}

void LocalSearch::SearchOrOpt2()
{
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		// improved = false;
		best_delta = 0;
		ry_idx = rx_idx;
		route_x = ads->routeADS[rx_idx];
		route_y = ads->routeADS[ry_idx];

		int counter = 0;

		while (counter < 100)
		{
			if (route_x->length <= 5 
					|| !curr_solution->GetIntraNeighbStatus(rx_idx, 2))
				break;

			improved = false;

#ifdef USEMMD
			intra_smd = ads->getSMD(rx_idx, 2);	
#endif

			start_search = 0;
			end_search  = route_x->length - 2;
			rx_ry_feasible = route_x->is_feasible;

			if (!use_mmd || intra_smd->status == UNVISITED)
			{
#ifdef USEMMD
				intra_smd->status = UNIMPROVED;
#endif
				for (pos_x = 1; pos_x < route_x->length - 2; pos_x++)
				{
					cl_x = curr_solution->route[rx_idx][pos_x];
					pos_z = pos_x + 1;
					cl_z = curr_solution->route[rx_idx][pos_z];
					prev_x = curr_solution->route[rx_idx][pos_x - 1];
					next_z = curr_solution->route[rx_idx][pos_z + 1];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_x];
						end_search = route_y->end_search[cl_x];
						// std::cout << "\n";
						// std::cout << start_search << " " << 0 << "\n";
						// std::cout << end_search << " " << route_y->length - 2 << "\n";
						// std::cout << "\n";
					}
#endif

					for (pos_y = start_search; pos_y <= end_search; pos_y++)
					{
						if (pos_y == pos_x || pos_y == pos_x + 1 || pos_x == pos_y + 1) 
							continue;
						cl_y = curr_solution->route[rx_idx][pos_y];
						pos_w = pos_y + 1;
						if (pos_x == pos_w || pos_y == pos_x)
							continue;
						cl_w = curr_solution->route[rx_idx][pos_w];
						// prev_y = curr_solution->route[rx_idx][pos_y - 1];
						EvalOrOpt();
					}
				}
			}
#ifdef USEMMD		
			if (intra_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = intra_smd->cost;
				best_pos_x = intra_smd->r1_pos1;
				best_pos_y = intra_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
			}
#endif

			if (improved)
			{
				stats->neighb_improv[8]++;
				//curr_solution->SetInterNeighbStatusAllTrue(rx_idx);
				//curr_solution->SetIntraNeighbStatusAllTrue(rx_idx);
				// std::cout << best_delta << "\n";
				// std::cout << *curr_solution << "\n";
				// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";
				// std::cout << "A\n";
				ExecRelocate2();
			}
			else
			{
				counter = 100;
				break;
		//		curr_solution->SetIntraNeighbStatus(rx_idx, 2, false);			
			}

			counter++;
		}
	}
}

void LocalSearch::SearchOrOpt3()
{
	improved = false;
	best_delta = 0;
	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		// improved = false;
		best_delta = 0;
		ry_idx = rx_idx;
		route_x = ads->routeADS[rx_idx];
		route_y = ads->routeADS[ry_idx];

		int counter = 0;

		while (counter < 100)
		{
			if (route_x->length <= 6
					|| !curr_solution->GetIntraNeighbStatus(rx_idx, 3))
				break;

			improved = false;

#ifdef USEMMD
			intra_smd = ads->getSMD(rx_idx, 3);	
#endif

			start_search = 0;
			end_search  = route_x->length - 2;
			rx_ry_feasible = route_x->is_feasible;

			if (!use_mmd || intra_smd->status == UNVISITED)
			{
#ifdef USEMMD
				intra_smd->status = UNIMPROVED;
#endif
				for (pos_x = 1; pos_x < (int) (route_x->length - 3); pos_x++)
				{
					cl_x = curr_solution->route[rx_idx][pos_x];
					pos_z = pos_x + 2;
					cl_z = curr_solution->route[rx_idx][pos_z];
					prev_x = curr_solution->route[rx_idx][pos_x - 1];
					next_z = curr_solution->route[rx_idx][pos_z + 1];

#ifdef USEFSR
					if (rx_ry_feasible)
					{
						start_search = route_y->start_search[cl_x];
						end_search = route_y->end_search[cl_x];
						// std::cout << "\n";
						// std::cout << start_search << " " << 0 << "\n";
						// std::cout << end_search << " " << route_y->length - 2 << "\n";
						// std::cout << "\n";
					}
#endif

					for (pos_y = start_search; pos_y < pos_x - 1; pos_y++)
					{
						cl_y = curr_solution->route[rx_idx][pos_y];
						pos_w = pos_y + 1;
						cl_w = curr_solution->route[rx_idx][pos_w];
						EvalOrOpt();
					}
					for (pos_y = pos_z + 1; pos_y <= end_search; pos_y++)
					{
						cl_y = curr_solution->route[rx_idx][pos_y];
						pos_w = pos_y + 1;
						cl_w = curr_solution->route[rx_idx][pos_w];
						EvalOrOpt();
					}
				}
			}
#ifdef USEMMD		
			if (intra_smd->cost + EPSILON < best_delta)
			{
				improved = true;
				best_delta = intra_smd->cost;
				best_pos_x = intra_smd->r1_pos1;
				best_pos_y = intra_smd->r2_pos1;
				best_route_x = route_x;
				best_route_y = route_y;
			}
#endif

			if (improved)
			{
				stats->neighb_improv[9]++;
				// std::cout << best_delta << "\n";
				// std::cout << "\n\n\n\n\n\n\n\n";
				// std::cout << *curr_solution << "\n";
				// std::cout << best_pos_x << " " << best_pos_y << " " << best_route_x << " " << best_route_y << "\n";
				// std::cout << "A\n";
				ExecRelocate3();
				//curr_solution->SetInterNeighbStatusAllTrue(rx_idx);
				//curr_solution->SetIntraNeighbStatusAllTrue(rx_idx);

				// std::cout << *curr_solution << "\n\n\n\n\n\n";
			}

			else
			{
				counter = 100;
				break;
				// curr_solution->SetIntraNeighbStatus(rx_idx, 3, false);			
			}

			counter++;
		}
	}
}

inline void LocalSearch::ExecRelocate3()
{
	if (best_route_x != best_route_y)
	{
		for (int i = best_pos_x; i < best_route_x->length - 2; i++)
		{
			std::swap(curr_solution->route[best_route_x->cour][i + 1], curr_solution->route[best_route_x->cour][i + 2]);
			std::swap(curr_solution->route[best_route_x->cour][i], curr_solution->route[best_route_x->cour][i + 1]);
		}
		int cust1, cust2;
		if (bestMoveType == 0)
		{
			cust2 = curr_solution->route[best_route_x->cour].back();
			curr_solution->route[best_route_x->cour].pop_back();
			cust1 = curr_solution->route[best_route_x->cour].back();
			curr_solution->route[best_route_x->cour].pop_back();
		}
		else
		{
			cust1 = curr_solution->route[best_route_x->cour].back();
			curr_solution->route[best_route_x->cour].pop_back();
			cust2 = curr_solution->route[best_route_x->cour].back();
			curr_solution->route[best_route_x->cour].pop_back();
		}

		best_route_x->Update();

		curr_solution->route[best_route_y->cour].push_back(cust1);
		curr_solution->route[best_route_y->cour].push_back(cust2);
		for (int i = best_route_y->length - 3; i > best_pos_y; i--)
		{
			std::swap(curr_solution->route[best_route_y->cour][i], curr_solution->route[best_route_y->cour][i + 1]);
			std::swap(curr_solution->route[best_route_y->cour][i + 1], curr_solution->route[best_route_y->cour][i + 2]);
		}
		best_route_y->Update();
	}

	// 0 1 2 3 4 5 6 7 8 9 0

	else
	{
		if (best_pos_x < best_pos_y)
		{
			std::rotate(curr_solution->route[best_route_x->cour].begin() + best_pos_x, curr_solution->route[best_route_x->cour].begin() + best_pos_x + 3, curr_solution->route[best_route_x->cour].begin() + best_pos_y + 1);
		}
		else
		{
			std::rotate(curr_solution->route[best_route_x->cour].begin() + best_pos_y + 1, curr_solution->route[best_route_x->cour].begin() + best_pos_x, curr_solution->route[best_route_x->cour].begin() + best_pos_x + 3);
		}
		best_route_x->Update();
	}

	UpdateClientRank();

	// std::cout << *curr_solution << "\n";
}


void LocalSearch::ResetThreeBestInsert()
{
	for (int i = 0; i < three_best_insert.size(); i++)
	{
		for (int j = 0; j < three_best_insert[i].size(); j++)
		{
			three_best_insert[i][j].best_cost[0] = 1e30;
			three_best_insert[i][j].best_cost[1] = 1e30;
			three_best_insert[i][j].best_cost[2] = 1e30;
		}
	}
}

void LocalSearch::SearchSwapStar()
{	
	improved = false;
	best_delta = 0;
	three_best_insert = std::vector<std::vector<ThreeBestInsertInfo>> (data->nb_clients + 1, std::vector<ThreeBestInsertInfo> (curr_solution->route.size()));
	// ResetThreeBestInsert();

	for (rx_idx = 0; rx_idx < ads->routeADS.size(); rx_idx++)
	{
		route_x = ads->routeADS[rx_idx];
		if (route_x->length <= 2)
			continue;

		for (ry_idx = rx_idx + 1; ry_idx < ads->routeADS.size(); ry_idx++)
		{
			if (rx_idx == ry_idx
				|| (!curr_solution->GetInterNeighbStatus(rx_idx, ry_idx, 6)))
				continue;
			route_y = ads->routeADS[ry_idx];

			rx_ry_improved  = false;

			ComputeThreeBestInsertions();			
			std::swap(rx_idx, ry_idx);
			std::swap(route_x, route_y);
			ComputeThreeBestInsertions();			
			std::swap(rx_idx, ry_idx);
			std::swap(route_x, route_y);

			for (pos_x = 1; pos_x <= route_x->length - 2; pos_x++)
			{
				cl_x  = curr_solution->route[rx_idx][pos_x];
				// TODO: switch prev_x, prev_y, ... with clients w and z accordingly
				int prev_x  = curr_solution->route[rx_idx][pos_x - 1];
				int next_x  = curr_solution->route[rx_idx][pos_x + 1];
				double removal_cost_x = - data->time_cost[prev_x][cl_x] - data->time_cost[cl_x][next_x] + data->time_cost[prev_x][next_x];

				for (pos_y = 1; pos_y <= route_y->length - 2; pos_y++)
				{		
					cl_y  = curr_solution->route[ry_idx][pos_y];
					int prev_y  = curr_solution->route[ry_idx][pos_y - 1];
					int next_y  = curr_solution->route[ry_idx][pos_y + 1];

					double delta_cap_viol_x = ComputeCapViolation(route_x->load.back() - data->client[cl_x].demand + data->client[cl_y].demand)
						- route_x->cap_violation;
					double delta_cap_viol_y = ComputeCapViolation(route_y->load.back() - data->client[cl_y].demand + data->client[cl_x].demand)
						- route_y->cap_violation;

					double removal_cost_y = - data->time_cost[prev_y][cl_y] - data->time_cost[cl_y][next_y] + data->time_cost[prev_y][next_y];

					ThreeBestInsertInfo* bti = &three_best_insert[cl_x][ry_idx];

					if (bti->prev[0] != cl_y && bti->next[0] != cl_y)
					{
						swapstar_ins_pos_x = bti->best_position[0];	
						swapstar_ins_cost_x = bti->best_cost[0];
					}
					else if (bti->prev[1] != cl_y && bti->next[1] != cl_y) 
					{
						swapstar_ins_pos_x = bti->best_position[1];	
						swapstar_ins_cost_x = bti->best_cost[1];
					}
					else 
					{
						swapstar_ins_pos_x = bti->best_position[2];	
						swapstar_ins_cost_x = bti->best_cost[2];
					}

					bti = &three_best_insert[cl_y][rx_idx];

					if (bti->prev[0] != cl_x && bti->next[0] != cl_x)
					{
						swapstar_ins_pos_y = bti->best_position[0];	
						swapstar_ins_cost_y = bti->best_cost[0];
					}
					else if (bti->prev[1] != cl_x && bti->next[1] != cl_x) 
					{
						swapstar_ins_pos_y = bti->best_position[1];	
						swapstar_ins_cost_y = bti->best_cost[1];
					}
					else 
					{
						swapstar_ins_pos_y = bti->best_position[2];	
						swapstar_ins_cost_y = bti->best_cost[2];
					}

					curr_delta = removal_cost_x + removal_cost_y + swapstar_ins_cost_x + swapstar_ins_cost_y 
						+ data->penalty_load * (delta_cap_viol_x + delta_cap_viol_y); 
					// std::cout << swapstar_ins_cost_x << " " << swapstar_ins_cost_y << "\n";
					// if (data->penalty_load * (delta_cap_viol_x + delta_cap_viol_y) < 0)
					// {
					// }
					//std::cout << route_x->cap_violation << "\n";
					//std::cout << route_y->cap_violation << "\n";

					rx_ry_improved |= (curr_delta + EPSILON < 0);
					if (curr_delta < best_delta)
					{
						improved = true;
						best_delta = curr_delta;
						best_pos_x = pos_x;
						best_pos_y = pos_y;
						best_swapstar_ins_pos_x = swapstar_ins_pos_x;
						best_swapstar_ins_pos_y = swapstar_ins_pos_y;
						best_route_x = route_x;
						best_route_y = route_y;
					}

				}
			}

		}
	}

	if (improved)
	{
		// std::cout <<"\n" << best_delta << "\n";

		stats->neighb_improv[11]++;
		ExecSwapStar();
		curr_solution->SetInterNeighbStatusAllTrue(best_route_x->cour, best_route_y->cour);
	}

	// 	std::cout <<"\n" << best_route_x->cour  << "\n";
	// 	std::cout <<"\n" << best_route_y->cour  << "\n";
	// 	std::cout <<"\n" << best_pos_x << "\n";
	// 	std::cout <<"\n" << best_pos_y << "\n";
	// 	std::cout <<"\n" << best_swapstar_ins_pos_x << "\n";
	// 	std::cout <<"\n" << best_swapstar_ins_pos_y << "\n";
	// 	int a = curr_solution->route[best_route_x->cour][best_pos_x];
	// 	int next_a = curr_solution->route[best_route_x->cour][best_pos_x + 1];
	// 	int prev_a = curr_solution->route[best_route_x->cour][best_pos_x - 1];
	// 
	// 	int b = curr_solution->route[best_route_y->cour][best_pos_y];
	// 	int next_b = curr_solution->route[best_route_y->cour][best_pos_y + 1];
	// 	int prev_b = curr_solution->route[best_route_y->cour][best_pos_y - 1];
	// 
	// 	int c = curr_solution->route[best_route_x->cour][best_swapstar_ins_pos_y];
	// 	int next_c = curr_solution->route[best_route_x->cour][best_swapstar_ins_pos_y + 1];
	// 
	// 	int d = curr_solution->route[best_route_y->cour][best_swapstar_ins_pos_x];
	// 	int next_d = curr_solution->route[best_route_y->cour][best_swapstar_ins_pos_x + 1];
	// 	int ra = best_route_x->cour;
	// 	int rb = best_route_y->cour;
	// 	std::vector<std::vector<double>> t = data->time_cost;
	// 
	// 	std::cout << - t[prev_a][a] - t[a][next_a] + t[prev_a][next_a]
	// 							- t[prev_b][b] - t[b][next_b] + t[prev_b][next_b]
	// 							+ t[c][b] + t[b][next_c] - t[c][next_c]
	// 							+ t[d][a] + t[a][next_d] - t[d][next_d] << "\n";
}

void LocalSearch::ComputeThreeBestInsertions()
{
	for (pos_x = 1; pos_x <= route_x->length - 2; pos_x++)
	{
		cl_x  = curr_solution->route[rx_idx][pos_x];
		for (pos_y = 0; pos_y <= route_y->length - 2; pos_y++)
		{		
			cl_y = curr_solution->route[ry_idx][pos_y];
			int next_y = curr_solution->route[ry_idx][pos_y + 1];
			double cost = - data->time_cost[cl_y][next_y] + data->time_cost[cl_y][cl_x] + data->time_cost[cl_x][next_y];
			three_best_insert[cl_x][ry_idx].CompareAndUpdateInsertion(pos_y, cost, cl_y, next_y);
		}
	}
}

void LocalSearch::ExecSwapStar()
{
	int best_cx = curr_solution->route[best_route_x->cour][best_pos_x];
	int best_cy = curr_solution->route[best_route_y->cour][best_pos_y];

	curr_solution->route[best_route_x->cour].erase(curr_solution->route[best_route_x->cour].begin() + best_pos_x);	
	curr_solution->route[best_route_y->cour].erase(curr_solution->route[best_route_y->cour].begin() + best_pos_y);	

	curr_solution->route[best_route_x->cour].insert(curr_solution->route[best_route_x->cour].begin() - (best_swapstar_ins_pos_y >= best_pos_x) + best_swapstar_ins_pos_y + 1, best_cy);	
	curr_solution->route[best_route_y->cour].insert(curr_solution->route[best_route_y->cour].begin() - (best_swapstar_ins_pos_x >= best_pos_y) + best_swapstar_ins_pos_x + 1, best_cx);	

	best_route_x->Update();
	best_route_y->Update();
}

/** perturbation procedures */

// void LocalSearch::multipleSwap()
// {
// 	std::vector<bool> modifiedRoutes(curr_solution->route.size(), false);
// 	int k = rand() % 2 + 1;
// 	for (int a = 0; a < k; a++)
// 	{
// 		double loadR1, loadR2;
// 		int r1, r2, i, j;
// 		do
// 		{
// 			loadR1 = loadR2 = 0;
// 			do
// 			{
// 				r1 = rand() % curr_solution->route.size() ;
// 				r2 = rand() % curr_solution->route.size() ;
// 			} while (curr_solution->route[r1].size() == 2 || curr_solution->route[r2].size() == 2 || r1 == r2);
// 			i = rand() % (curr_solution->route[r1].size() - 2) + 1;
// 			j = rand() % (curr_solution->route[r2].size() - 2) + 1;
// 			for (int i = 0; i < curr_solution->route[r1].size(); i++)
// 			{
// 				loadR1 += data->client[curr_solution->route[r1].demand[i]];
// 			}
// 			for (int i = 0; i < curr_solution->route[r2].size(); i++)
// 			{
// 				loadR2 += data->client[curr_solution->route[r2].demand[i]];
// 			}
// 			// loadR1 += (-data->client[curr_solution->route[r1].demand[i]] + data->client[curr_solution->route[r2].demand[j]]);
// 			// loadR2 -= (-data->client[curr_solution->route[r1].demand[i]] + data->client[curr_solution->route[r2].demand[j]]);
// 		} while (loadR1 > data->capacity || loadR2 > data->capacity);
// 
// 		std::swap(curr_solution->route[r1][i], curr_solution->route[r2][j]);
// 
// 		// SWITCH TO VECTOR TO UPDATE LATER
// 		modifiedRoutes[r1] = true;
// 		modifiedRoutes[r2] = true;
// 	}
// 
// 	for (int r = 0; r < modifiedRoutes.size(); r++)
// 	{
// 		if (modifiedRoutes[r])
// 			routeData[r].Update();
// 	}
// 
// 	UpdateClientRank();
// 	// ComputePenalizedCost(curr_solution);
// }
// 
// void LocalSearch::multipleShift()
// {
// 	bool perturbed = true;
// 	int num_shifts = rand() % 2 + 1;
// 
// 	std::vector<bool> modifiedRoutes(curr_solution->route.size(), false);
// 
// 	for (int shift_counter = 0; shift_counter < num_shifts; shift_counter++)
// 	{
// 		int r1 = rand() % ads->routeADS.size();
// 		while (routeData[r1].nbNodes == 2)
// 		{
// 			r1 = rand() % ads->routeADS.size();
// 		}
// 		int r2 = rand() % ads->routeADS.size();
// 		while (r2 == r1 || routeData[r2].nbNodes == 2)
// 			r2 = rand() % ads->routeADS.size();
// 
// 		int r1_size = routeData[r1].nbNodes;
// 		int r2_size = routeData[r2].nbNodes;
// 
// 		int k_index = 1 + rand() % (routeData[r1].nbNodes - 2);
// 		int l_index = 1 + rand() % (routeData[r2].nbNodes - 2);
// 		int k = routeData[r1].path[k_index];
// 		int l = routeData[r2].path[l_index];
// 
// 		// s->Print();
// 		// cout << r1 << " " << r2 << "\n";
// 		// cout << k_index << " " << l_index << "\n";
// 		// cout << k << " " << l << "\n";
// 
// 		{
// 			perturbed = true;
// 
// 			int r1_index = 1 + rand() % (routeData[r1].path.size() - 2);
// 			int r2_index = 1 + rand() % (routeData[r2].path.size() - 2);
// 
// 			routeData[r1].path.erase(routeData[r1].path.begin() + k_index);
// 			routeData[r2].path.erase(routeData[r2].path.begin() + l_index);
// 
// 			routeData[r1].path.insert(routeData[r1].path.begin() + r1_index, l);
// 			routeData[r2].path.insert(routeData[r2].path.begin() + r2_index, k);
// 
// 			modifiedRoutes[r1] = true;
// 			modifiedRoutes[r2] = true;
// 
// 			// s->cost = s->CalculateCost();
// 		}
// 	}
// 
// 	for (int r = 0; r < modifiedRoutes.size(); r++)
// 	{
// 		if (modifiedRoutes[r])
// 			routeData[r].Update();
// 	}
// 
// 	UpdateClientRank();
// 	// ComputePenalizedCost(curr_solution);
// }
