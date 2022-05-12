#include "IlsRvnd.h"

AdaptativeInfo::AdaptativeInfo(
		size_t max_nb_solutions, 
		double decrease_factor
) : max_nb_solutions(max_nb_solutions), 
	decrease_factor(decrease_factor) {}

void AdaptativeInfo::Update(Solution *solution)
{
	UpdatePrevSolutions(solution);
	UpdateBestSolutionIdx();
	UpdateQualityFactor();
	UpdateThreshold(solution);

	if (solution->cost < threshold)
	{
		reference_sol_idx = last_sol_idx;
	}

	return;
}

void AdaptativeInfo::UpdatePrevSolutions(Solution *solution)
{
	if (previous_solutions.size() < max_nb_solutions)
	{
		previous_solutions.push_back(*solution);
	} 
	else 
	{
		previous_solutions[last_sol_idx] = *solution;
		last_sol_idx++;
		if (last_sol_idx == max_nb_solutions)
		{
			last_sol_idx = 0;
		}
	}

	return;
}

void AdaptativeInfo::UpdateBestSolutionIdx()
{
	double best_cost = std::numeric_limits<double>::infinity();

	for (int i = 0; i < previous_solutions.size(); i++)
	{
		if (previous_solutions[i].cost < best_cost)
		{
			best_cost = previous_solutions[i].cost;
			best_sol_idx = i;
		}
	}

	return;
}

void AdaptativeInfo::UpdateQualityFactor()
{
	quality_factor *= decrease_factor;
	// cout << quality_factor << endl;

	return;
}

void AdaptativeInfo::UpdateThreshold(Solution *solution)
{
	int nb_solutions = previous_solutions.size();
	average_cost = average_cost * (1.0 - (1.0 / nb_solutions))
								+ ((double) solution->cost / nb_solutions);

	double best_cost = previous_solutions[best_sol_idx].cost;

	threshold = best_cost + 
							quality_factor * (average_cost - best_cost);

	// cout << threshold << endl;

	return;
}

IlsRvnd::IlsRvnd(Data *data) : data(data)
{
}

IlsRvnd::IlsRvnd(Data *data, Pool *my_pool, Solution *s0) : data(data), my_pool(my_pool)
{
	init_solution = s0;
	skip_construction = true;
	my_sp = NULL;
}

IlsRvnd::IlsRvnd(Data *data, Pool *my_pool) : data(data), my_pool(my_pool)
{
	my_sp = new SPSolverCplex(data, my_pool, &best);
}

void IlsRvnd::LogCurrentIterationStats(IlsStats *stats)
{
	// std::vector<std::string> strVec =
	// 	{"Shift(1,0)",
	// 	 "Shift(2,0)",
	// 	 "Swap(1,1)",
	// 	 "Swap(2,1)",
	// 	 "Swap(2,2)",
	// 	 "Cross",
	// 	 "Exchange",
	// 	 "Reinsertion",
	// 	 "Or-opt-2",
	// 	 "Or-opt-3",
	// 	 "2-OPT"};

	// std::cout << std::endl
	// 		  << "************************************************************************************************************************************************" << std::endl
	// 		  << std::endl;
	// std::cout << " Iteration: " << iteration << std::endl;
	// std::cout << " Improvements: " << ils_improv << std::endl;
	// std::cout << " Initial sol cost: " << init_sol_cost << std::endl;
	// std::cout << " Best solution cost: " << best.cost << std::endl;
	// std::cout << " Last Iteration improvement: " << stats->last_improv_iter << std::endl
	// 		  << std::endl;
	// std::cout << " Neighborhood:    ";
	// for (int i = 0; i < strVec.size(); i++)
	// {
	// 	std::cout << std::setw((int)strVec[i].size() + 3) << std::right << strVec[i];
	// }
	// std::cout << std::endl;
	// std::cout << " Nb improvements: ";
	// for (int i = 0; i < strVec.size(); i++)
	// {
	// 	std::cout << std::setw((int)strVec[i].size() + 3) << std::right << std::to_string(stats->neighb_improv[i]);
	// }

	// std::cout << std::endl
	// 		  << std::endl;
	// std::cout << "************************************************************************************************************************************************" << std::endl;

	std::vector<std::string> strVec =
		{"Sh1",
		 "Sh2",
		 "Sw11",
		 "Sw21",
		 "Sw22",
		 "C",
		 "Ex",
		 "Or1",
		 "Or2",
		 "Or3",
		 "2opt",
		 "Sw*"};

	std::vector<std::string> pertVec =
		{"P_Sw11",
		"P_Sh11"};

	std::cout << "\n";
	//std::cout << "last improv: " << stats->last_improv_iter << " ";
	//std::cout << "iter " << iteration << ", ";
	std::cout << ils_improv << " improvs, ";
	std::cout << init_sol_cost << " -> ";
	std::cout << best.cost << " ";

	for (int i = 0; i < strVec.size(); i++)
	{
		std::cout << strVec[i] << "-" << std::to_string(stats->neighb_improv[i]) << " ";
	}
	for (int i = 0; i < pertVec.size(); i++)
	{
		std::cout << pertVec[i] << "-" << std::to_string(stats->perturb_improv[i]) << " ";
	}

	std::cout << std::endl;
}

void IlsRvnd::LogFinalStats(IlsStats *stats)
{
	puts("******************************************************************************************************************************************************");
	std::cout << std::setw(19) << std::right << "Vizinhanca"
			  << std::setw(19) << std::right << "Tempo Hash"
			  << std::setw(19) << std::right << "Consultas Hash"
			  << std::setw(19) << std::right << "Templo Explor."
			  << std::setw(19) << std::right << "Evitadas";
	std::cout << std::endl;

	std::cout << std::setw(19) << std::right << "Shift1"
			  << std::setw(19) << std::right << stats->inter_mmd_search[0].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_search[0].count
			  << std::setw(19) << std::right << stats->inter_neighb_explore[0].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_effective[0];
	std::cout << std::endl;

	std::cout << std::setw(19) << std::right << "Shift2"
			  << std::setw(19) << std::right << stats->inter_mmd_search[4].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_search[4].count
			  << std::setw(19) << std::right << stats->inter_neighb_explore[4].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_effective[4];
	std::cout << std::endl;

	std::cout << std::setw(19) << std::right << "Swap11"
			  << std::setw(19) << std::right << stats->inter_mmd_search[1].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_search[1].count
			  << std::setw(19) << std::right << stats->inter_neighb_explore[1].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_effective[1];
	std::cout << std::endl;

	std::cout << std::setw(19) << std::right << "Swap22"
			  << std::setw(19) << std::right << stats->inter_mmd_search[2].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_search[2].count
			  << std::setw(19) << std::right << stats->inter_neighb_explore[2].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_effective[2];
	std::cout << std::endl;

	std::cout << std::setw(19) << std::right << "Swap21"
			  << std::setw(19) << std::right << stats->inter_mmd_search[5].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_search[5].count
			  << std::setw(19) << std::right << stats->inter_neighb_explore[5].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_effective[5];
	std::cout << std::endl;

	std::cout << std::setw(19) << std::right << "Cross"
			  << std::setw(19) << std::right << stats->inter_mmd_search[3].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_search[3].count
			  << std::setw(19) << std::right << stats->inter_neighb_explore[3].total_time
			  << std::setw(19) << std::right << stats->inter_mmd_effective[3];
	std::cout << std::endl;

	puts("******************************************************************************************************************************************************");
}

int IlsRvnd::EstimateNumberOfVehicles()
{
	int v = 1;
	std::vector<int> CL(data->nb_clients);
	
	for (int i = 1; i <= data->nb_clients; i++)
		CL[i - 1] = i;

	std::swap(CL[rand() % CL.size()], CL.back());
	int k = CL.back();
	CL.pop_back();

	std::vector<std::vector<int>> s(1, std::vector<int>{0, k, 0});
	std::vector<double> load(1, data->client[k].demand);
	auto best_node = CL.begin();
	int best_pos = 0;
	bool force_feasible = true;

	while (!CL.empty())
	{
		bool found_client = false;
		double best_insertion = std::numeric_limits<double>::infinity();
		int best_route;
		for (auto p = CL.begin(); p != CL.end(); ++p)
		{
			int r = s.size() - 1;
			for (int pos = 0; pos < s[r].size() - 1; pos++)
			{
				double insertion_cost = data->time_cost[s[r][pos]][*p] + data->time_cost[*p][s[r][pos + 1]] - data->time_cost[s[r][pos]][s[r][pos + 1]];

				if (insertion_cost < best_insertion)
				{
					double total_demand = load[r] + data->client[*p].demand;

#ifdef ROBUST_DEMAND
					std::vector<double> sorted_deviations(s[r].size() - 1);

					for (int i = 1; i < s[r].size() - 1; i++)
					{
						sorted_deviations[i - 1] = data->client[s[r][i]].demand_deviation;
						// total_demand += data->client[s[r][i]].demand;
					}
					sorted_deviations.back() = data->client[s[r][pos]].demand_deviation;

					if (data->budget_demand > 0)
					{
						if (s[r].size() - 2 >= data->budget_demand)
							std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin(), sorted_deviations.end());
						for (int i = 0; i < sorted_deviations.size() && i < data->budget_demand; i++)
							total_demand += sorted_deviations[i];
					}
#endif

					if (total_demand > data->capacity && force_feasible)
						continue;

					best_insertion = insertion_cost;
					best_node = p;
					best_pos = pos + 1;
					found_client = true;
					best_route = r;
				}
			}
		}

		if (found_client)
		{
			s[best_route].insert(s[best_route].begin() + best_pos, *best_node);
			load[best_route] += data->client[(*best_node)].demand;
			CL.erase(best_node);
		}
		else
		{
			s.push_back({0, 0});
			load.push_back(0);
		}
	}

	for (int i = 0; i < s.size(); i++)
	{
		double totdem = 0;
		for (int j = 0; j < s[i].size(); j++)
		{
		// 	std::cout << s[i][j] << " ";

			if (j > 0 && j < s[i].size() - 1)
				totdem += data->client[s[i][j]].demand;
		}
		//if (std::abs(load[i] - totdem) > 0.01)
		//{
		//	std::cout << "PERIGO\n";
		//	exit(0);
		//}
		//std::cout << load[i] << " " << totdem << " " << data->capacity;
		//std::cout << std::endl;
	}

	return (int) s.size();
}

void IlsRvnd::Run(int maxIter, int max_iter_ils, double tol)
{
	tolerance = tol;
	int temp_pool_counter = 0;

	//if(!skip_construction)
	{
		best = Solution(data);
		best.cost = std::numeric_limits<double>::max();
	}

	ADS *ads = new ADS(data);
	Solution *curr_solution = new Solution(data);
	Solution *best_iter_solution = new Solution;
	Solution *ref_solution = new Solution;
	Perturbation *perturb = new Perturbation(curr_solution, ads);
	LocalSearch *local_search = new LocalSearch(curr_solution, ads);
	IlsStats *stats = new IlsStats;

	int vehi_estimate;
	if (!skip_construction)
	{
		if (data->variant == "rcvrp")
			vehi_estimate = EstimateNumberOfVehicles();
		else
			vehi_estimate = data->nb_vehicles;

		std::cout << "Nb vehicles: " << vehi_estimate << "\n";
		data->nb_vehicles = vehi_estimate;
		data->max_iter_ils = data->nb_clients + 5 * vehi_estimate;
		max_iter_ils = data->nb_clients + 5 * vehi_estimate;

		if (data->nb_clients <= 150 && (data->nb_clients / vehi_estimate) < 11)
		{
			data->tolerance = 0.05;
		}
		else if (data->nb_clients > 150)
		{
			if (data->nb_clients / vehi_estimate < 11)
				data->tolerance = 0.005;
			else
				data->tolerance = 0.05;
			// max_iter_ils = 1500;
		}

		data->tolerance = 1;
	}

	int conta = 0;
	iteration = 0;
	//for (iteration = 0; iteration < maxIter; iteration++)
	//while (iteration++ < maxIter)
	// while (iteration++ < maxIter)
	double best_ils_cost = std::numeric_limits<double>::infinity();
	int when_to_clean = 2;
	int sp_step = 2;
	int sp_call_count = 0;
	
	while (iteration++ < maxIter)
	{

		// if (!skip_construction)
		// 	std::cout << std::endl << " Iter " << iteration << std::endl;
		auto start = std::chrono::system_clock::now();

		ils_improv = 0;
		ads->inter_mov_memo.clear();
		ads->intra_mov_memo.clear();
		local_search->stats = stats;
		// std::cout << "Iteration " << iteration << "\n";
		// data->nb_vehicles = 49;

		if (!skip_construction)
		{
			*curr_solution = GenerateInitialSolution(data, maxIter, vehi_estimate);
			vehi_estimate = curr_solution->route.size();
		}
		else
		{
			*curr_solution = *init_solution;
		}

		// std::cout << "Número de veículos: " << curr_solution->nb_vehicles << "\n";

		ads->SetSolution(curr_solution);
		ads->SetPenalizedCost();
		init_sol_cost = curr_solution->cost;

		*best_iter_solution = *curr_solution;
		*ref_solution = *curr_solution;

		cost_history = std::deque<std::tuple<double,double>>();


		// ILS procedure starts here
		AdaptativeInfo *adap = new AdaptativeInfo(30, 0.9);
		int iter_ils = 0;
		int nb_it_ils = 0;

		while (iter_ils < max_iter_ils)
		{
			iter_id++;
			conta++;
			ads->SetSolution(curr_solution);
			local_search->Run();
			ads->SetPenalizedCost();

#ifdef SP_STATS
			curr_solution->route_stats.assign(curr_solution->route.size(), 
					RouteStats(ils_improv, (skip_construction ? -1 : iteration), 
							false, curr_solution->cost, best_of_all_cost));
#endif

#ifdef USEMSD
			if (ads->FindAndPush(iteration))
			{
				// std::cout << "Found solution\n";
				break;
			}
#endif

			if (my_pool && (curr_solution->cost <= (1 + data->tolerance) * best.cost))
			{
					my_pool->Update(curr_solution, false);
			}

			cost_history.push_back({curr_solution->cost, best_iter_solution->cost});
			auto end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;

			if (curr_solution->cost + EPSILON < best_iter_solution->cost)
			{
				ils_improv++;
				*best_iter_solution = *curr_solution;
				iter_ils = 0;

				//std::cout << "(" << best_iter_solution->cost << ", " << std::flush;
				//std::cout << std::fixed << elapsed_seconds.count() << std::setprecision(2) << ") " << std::flush;
				// printf("(%g, %.2f) ", best_iter_solution->cost, 100.0 * (best_iter_solution->cost - best_of_all_cost) / best_of_all_cost);
				//printf("(%g, %.2f) ", best_iter_solution->cost, elapsed_seconds.count());
				// std::cout << std::flush;
				
				if (perturb->last_perturb != -1)
					stats->perturb_improv[perturb->last_perturb]++;
			}

			// adap->Update(curr_solution);

			// std::cout << data->tolerance << "\n";
			// adap->UpdatePerturbDegree(curr_solution, ref_solution);

			// *curr_solution = adap->previous_solutions[adap->reference_sol_idx];
			*curr_solution = *best_iter_solution;

			ads->SetSolution(curr_solution);
			//perturb->Run(adap->perturb_degree);
			perturb->Run(0);
			// std::cout << adap->perturb_degree << "\n";
			iter_ils++;

#ifdef DIMACS
//			if (best_iter_solution->cost < best.cost)
		{
			if (global_solution)
			{
				global_solution->partial_solution[*clust_idx] = *best_iter_solution;
				global_solution->Print();
			}
			else
				std::cout << *best_iter_solution << "\n";
		}
#endif
		}

		ils_gap = (double) (best_iter_solution->cost - best_ils_cost) / best_ils_cost * 100;
		sp_gap = (double) (best_iter_solution->cost - best.cost) / best.cost * 100;

		if (best_iter_solution->cost < best_ils_cost)
		{
			best_ils_cost = best_iter_solution->cost;
		}

		if (best_iter_solution->cost < best.cost)
		{
			best = *best_iter_solution;
			best_of_all_cost = best.cost;
			LogCurrentIterationStats(stats);
		}
		else
		{
			std::cout << "." << std::flush;
		}
		stats->ResetNeighbImprov();

		// if (!skip_construction && iteration > 1)
		// {
		// 	std::cout << "spG: " << sp_gap << "% "
		// 		<< " ilsG: " << 	 ils_gap << "% "
		// 		<< std::endl;
		// }

		if (!skip_construction)
		{
			history_counter = 0;

			std::ofstream history_file;
			history_file.open("plotdata/cost_iter" + std::to_string(iteration) + ".dat", std::ios::out);

			int i = 0;
			for (auto it = cost_history.begin(); it != cost_history.end(); ++it, i++)
			{
				history_file << i << " " << std::get<0>(*it) << " "
							 << std::get<1>(*it) << " " << best_iter_solution->cost << " "
							 << best.cost << "\n";
			}
			history_file.close();
		}

		// if (true && !skip_construction && (data->nb_clients >= 150) && my_sp)
		if (true && !skip_construction && (data->nb_clients >= 150) && my_sp)
		{
			bool sp_improved = true;
			RouteStats stats(ils_improv, (skip_construction ? -1 : iteration), true, 
					curr_solution->cost, best_of_all_cost);
			my_pool->Update(best_iter_solution, true);
			my_sp->SetSolution(&best);

			while (sp_improved)
			{
				//std::cout << "iter : " << iteration << "\n";
				std::cout << std::endl << "Perm pool size: " << my_pool->GetNbPermPath() << std::endl;
				std::cout << "Temp pool size: " << my_pool->GetNbNonPermPath() << std::endl;
				std::cout << "Max gap: " << my_pool->gap_tol << std::endl;

//				my_pool->GapClear(best_of_all_cost);
				my_sp->Run();

				if (my_sp->getBestSolution().cost < best.cost)
				{
					sp_improved = true;
					best = my_sp->getBestSolution();
					best_of_all_cost = best.cost;
					my_sp->SetSolution(&best);
					my_pool->Update(&best, true);
					iteration++;
					temp_pool_counter++;
				}
				else
				{
					sp_improved = false;
				}
			}

			temp_pool_counter++;

			if (temp_pool_counter >= when_to_clean)
			{
				std::cout << "\nCleaning temporary pool...\n";
				my_pool->RemoveTempRoutes();
				temp_pool_counter = 0;
			}
		}
		
		delete adap;
	}
	std::cout << std::endl;

	// LogFinalStats(stats);

	// Set partitioning
	if (true && !skip_construction && (data->nb_clients < 150) && my_sp)
	{
		double prev_cost;
		do
		{
			// my_sp->bestSolution = *best;
			std::cout << "Perm pool size: " << my_pool->GetNbPermPath() << "\n";
			std::cout << "Temp pool size: " << my_pool->GetNbNonPermPath() << "\n";
			prev_cost = best.cost;
			my_pool->Update(&best, true);
			my_sp->SetSolution(&best);
			my_sp->Run();
			best = my_sp->getBestSolution();
			best_of_all_cost = best.cost;
		} while (best.cost + EPSILON < prev_cost);
	}

	delete stats;
	delete ads;
	delete curr_solution;
	delete best_iter_solution;
	delete ref_solution;
	delete perturb;
	delete local_search;
	//	delete my_sp;
	//	delete best;

	//std::cout << "\ncontador:" << ads->contador1 << "\n";
	//std::cout << conta << "\n\n\n\n";
}
