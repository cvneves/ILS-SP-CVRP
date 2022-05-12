#include "Construction.h"

bool ParallelInsertion(Data *data, Solution *solution)
{
    /** Initializing the candidate list CL with all clients */

    std::vector<int> CL(data->nb_clients);

    for (int i = 1; i <= data->nb_clients; i++)
        CL[i - 1] = i;

    std::random_shuffle(CL.begin(), CL.end());

    solution->route = std::vector<std::vector<int>>(solution->nb_vehicles);
    solution->route_cost = std::vector<double>(solution->nb_vehicles);
    solution->route_load = std::vector<double>(solution->nb_vehicles);
    std::vector<double> load = std::vector<double>(solution->nb_vehicles, 0.0);

    /** Assigning a single client from CL to each route randomly */
    for (int r = 0; r < solution->nb_vehicles; r++)
    {
        if (CL.empty())
        {
            solution->route[r] = {0, 0};
            continue;
        }
				if (r < solution->nb_vehicles - 1)
				{
        	solution->route[r] = {0, CL.back(), 0};
        	load[r] = data->client[CL.back()].demand;
        	CL.pop_back();
				}
				else 
				{
        	solution->route[r] = {0, 0};
				}
    }

    int cheapestInsertion = rand() % 2;

    bool forceFeasible = true;

    double gamma = (rand() % 35);
    gamma = gamma / 20.0;

    /** Pick a client and insert it into the best position possible until CL is empty */
    while (!CL.empty())
    {
        bool foundClient = false;
        double bestInsertionCost = std::numeric_limits<double>::infinity();
        auto bestNode = CL.begin();
        int bestroute;
        int bestPosition;

        for (auto p = CL.begin(); p != CL.end(); ++p)
        {
            for (int r = 0; r < solution->nb_vehicles; r++)
            {
                for (int pos = 0; pos < solution->route[r].size() - 1; pos++)
                {
                    double insertionCost;

                    if (cheapestInsertion)
                    {
                        insertionCost = data->time_cost[solution->route[r][pos]][*p]
												 	+ data->time_cost[*p][solution->route[r][pos + 1]]
												 	- data->time_cost[solution->route[r][pos]][solution->route[r][pos + 1]] 
													- gamma * (data->time_cost[0][*p] + data->time_cost[*p][0]);
                    }
                    else
                    {
                        insertionCost = data->time_cost[solution->route[r][pos]][*p];
                    }

                    if (insertionCost < bestInsertionCost)
                    {
                        double totalDemand = load[r] + data->client[*p].demand;
                        std::vector<double> sorted_deviations(solution->route[r].size() - 1);
                        for (int i = 1; i < solution->route[r].size() - 1; i++)
                        {
                            sorted_deviations[i - 1] = data->client[solution->route[r][i]].demand_deviation;
                            // totalDemand += data->client[route[r][i]].demand;
                        }
                        sorted_deviations.back() = data->client[solution->route[r][pos]].demand_deviation;

                        if (data->budget_demand > 0)
                        {
                            if (solution->route[r].size() - 2 >= data->budget_demand)
                                std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin(), sorted_deviations.end());
                            for (int i = 0; i < sorted_deviations.size() && i < data->budget_demand; i++)
                                totalDemand += sorted_deviations[i];
                        }

                        if (totalDemand > data->capacity && forceFeasible)
                            continue;
                        foundClient = true;
                        bestInsertionCost = insertionCost;
                        bestNode = p;
                        bestroute = r;
                        bestPosition = pos + 1;
                    }
                }
            }
        }

        if (foundClient)
        {
            solution->route[bestroute].insert(solution->route[bestroute].begin() + bestPosition, *bestNode);
            load[bestroute] += data->client[(*bestNode)].demand;
            CL.erase(bestNode);
        }
        else
        {
#ifdef ROBUST_DEMAND
            forceFeasible = false;
#else
            break;
#endif
        }
    }

    if (CL.empty())
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool SequentialInsertion(Data *data, Solution *solution)
{
    /** Initializing the candidate list CL with all clients */

    std::vector<int> CL(data->nb_clients);

    for (int i = 1; i <= data->nb_clients; i++)
        CL[i - 1] = i;

    std::random_shuffle(CL.begin(), CL.end());

    solution->route = std::vector<std::vector<int>>(solution->nb_vehicles);
    solution->route_cost = std::vector<double>(solution->nb_vehicles);
    solution->route_load = std::vector<double>(solution->nb_vehicles);
    std::vector<double> load(solution->nb_vehicles, 0.0);

    /** Assigning a single client from CL to each solution->route randomly */
    for (int r = 0; r < solution->nb_vehicles; r++)
    {
        if (CL.empty())
        {
            solution->route[r] = {0, 0};
            continue;
        }
        solution->route[r] = {0, CL.back(), 0};
        load[r] = data->client[CL.back()].demand;
        CL.pop_back();
    }

    int cheapestInsertion = rand() % 2;

    bool forceFeasible = true;
    std::vector<int> routeAvailable(solution->nb_vehicles, 1);
    int nbFullRoutes = 0;

    double gamma = (rand() % 35);
    gamma = gamma / 20.0;

    /** Pick a client and insert it into the best position possible until CL is empty */
    int r = 0;
    while ((nbFullRoutes < solution->nb_vehicles || !forceFeasible) && !CL.empty())
    {
        bool foundClient = false;
        double bestInsertionCost = std::numeric_limits<double>::infinity();
        auto bestNode = CL.begin();
        int bestPosition;

        if (routeAvailable[r] == 1)
        {
            for (int pos = 0; pos < solution->route[r].size() - 1; pos++)
            {
                for (auto p = CL.begin(); p != CL.end(); ++p)
                {
                    double insertionCost;

                    if (cheapestInsertion)
                    {
                        insertionCost = data->time_cost[solution->route[r][pos]][*p] + data->time_cost[*p][solution->route[r][pos + 1]] - data->time_cost[solution->route[r][pos]][solution->route[r][pos + 1]] -
                                        gamma * (data->time_cost[0][*p] + data->time_cost[*p][0]);
                    }
                    else
                    {
                        insertionCost = data->time_cost[solution->route[r][pos]][*p];
                    }

                    if (insertionCost < bestInsertionCost)
                    {
                        double totalDemand = load[r] + data->client[*p].demand;
                        std::vector<double> sorted_deviations(solution->route[r].size() - 1);
                        for (int i = 1; i < solution->route[r].size() - 1; i++)
                        {
                            sorted_deviations[i - 1] = data->client[solution->route[r][i]].demand_deviation;
                            // totalDemand += data->client[solution->route[r][i]].demand;
                        }
                        sorted_deviations.back() = data->client[solution->route[r][pos]].demand_deviation;

#ifdef ROBUST_DEMAND
                        if (data->budget_demand > 0)
                        {
                            if (solution->route[r].size() - 2 >= data->budget_demand)
                                std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin(), sorted_deviations.end());
                            for (int i = 0; i < sorted_deviations.size() && i < data->budget_demand; i++)
                                totalDemand += sorted_deviations[i];
                        }
#endif

                        if (totalDemand > data->capacity && forceFeasible)
                            continue;

                        foundClient = true;
                        bestInsertionCost = insertionCost;
                        bestNode = p;
                        bestPosition = pos + 1;
                    }
                }
            }

            if (foundClient)
            {
                solution->route[r].insert(solution->route[r].begin() + bestPosition, *bestNode);
                load[r] += data->client[(*bestNode)].demand;
                CL.erase(bestNode);
            }
            else
            {
                routeAvailable[r] = 0;
                nbFullRoutes++;
#ifdef ROBUST_DEMAND
                if (nbFullRoutes == solution->nb_vehicles)
                {
                    forceFeasible = false;
                }
#endif
            }
        }
        if (r == solution->nb_vehicles - 1)
        {
            r = -1;
        }
        r++;
    }

    if (CL.empty())
    {
        return true;
    }
    else
    {
        return false;
    }
}

Solution GenerateInitialSolution(Data *data, int maxIter, int nb_vehicles)
{
    Solution solution(data);
    int consecutiveTrials = 0;
    bool foundSolution;
    solution.nb_vehicles = nb_vehicles;
		
    while (true)
    {
#ifndef ROBUST_DEMAND
				int constructionType = rand() % 2;
#else
				int constructionType = 0;
#endif
        if (constructionType)
        {
            foundSolution = SequentialInsertion(data, &solution);
        }
        else
        {
            foundSolution = ParallelInsertion(data, &solution);
        }

        if (foundSolution)
        {
					solution.InitializeNeighbStatus();
					return solution;
        }
        else
        {
            consecutiveTrials++;
            if (consecutiveTrials == maxIter && !data->limit_nb_vehicles)
            {
                solution.nb_vehicles++;
								std::cout << "Nb. of vehicles updated to: " << solution.nb_vehicles << "\n";
                consecutiveTrials = 0;

								best_nb_vehi = std::max(solution.nb_vehicles, best_nb_vehi);
            }
        }
    }
}
