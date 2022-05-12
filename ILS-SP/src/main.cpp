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
#include <iomanip>
#include <map>

#include "ArgParser.h"
#include "Data.h"
#include "Solution.h"
#include "LocalSearch.h"
#include "IlsRvnd.h"
#include "SP.h"
#include "Construction.h"

// #include <mlpack/core.hpp>
// #include <mlpack/methods/kmeans/kmeans.hpp>
// 
// using namespace mlpack;
// using namespace mlpack::kmeans;
// using namespace mlpack::util;

// ./rvrp.out -i instances/TW/100/C105.txt -var rvrptw -bt 0 -at 0 -bq 5 -aq 0.5 
// ./rvrp.out -i instances/CVRP/Vrp-Set-X/X/X-n106-k14.vrp -aq 0.25 -bq 5 -at 0.1 -bt 0 -var rcvrp

double IlsRvnd::best_of_all_cost = 999999999;
int IlsRvnd::iter_id = 0;

#ifndef TESTING
bool RunSavings(Solution &best_of_all, Data &data);

int main(int argc, char *argv[])
{
	ArgParser arg(argc, argv);

	Data data(arg.params);

	// if (false && data.nb_clients > 500 && data.dist_type != 2)
	// {

	// 	ArgParser arg(argc, argv);

	// 	Data data(arg.params);

	// 	if (data.nb_clients > 500)
	// 	{
	// 		auto start = std::chrono::system_clock::now();

	// 		arma::mat points(2, data.nb_clients);
	// 		for (int i = 1; i <= data.nb_clients; i++)
	// 		{
	// 			points(0, i - 1) = data.client[i].x;
	// 			points(1, i - 1) = data.client[i].y;
	// 		}
	// 		arma::Row<size_t> assignments;
	// 		arma::mat centroids;
	// 		size_t nb_clusters = data.nb_clients / 200;

	// 		KMeans<> k;
	// 		k.Cluster(points, nb_clusters, assignments, centroids);

	// 		arma::Row<size_t> super_assignments;
	// 		// arma::mat points_centroids(2, nb_clusters);
	// 		// for (int i = 0, j = 0; i < nb_clusters; i++, j+=2)
	// 		// {
	// 		// 	points_centroids(0, i) = centroids[j];
	// 		// 	points_centroids(1, i) = centroids[j+1];
	// 		// }
	// 		arma::mat super_centroids;
	// 		size_t nb_super_clusters = std::max(1, (int) (nb_clusters / 2.5));
	// 		KMeans<> super_k;
	// 		super_k.Cluster(centroids, nb_super_clusters, super_assignments, super_centroids);

	// 		//std::cout << "digraph G {\n";

	// 		//std::vector<std::string> colors = {"blue",
	// 		//	"chartreuse",
	// 		//	"cadetblue",
	// 		//	"chocolate",
	// 		//	"darkorchid",
	// 		//	"deeppink",
	// 		//	"coral",
	// 		//	"mediumpurple",
	// 		//	"red"};

	// 		std::map<int, int> clust_count;
	// 		std::vector<std::vector<int>> cluster(nb_clusters, std::vector<int>());
	// 		for (int i = 1; i <= data.nb_clients; i++)
	// 		{
	// 			//std::cout << "v" << i << " [shape = point, pos=\"" << data.client[i].x / 1000.0 << "," << data.client[i].y / 1000.0 << "!\"";

	// 			//std::cout << ", color=";
	// 			//std::cout << colors[(assignments[i - 1] - 1) % colors.size()];
	// 			cluster[assignments[i - 1]].push_back(i);
	// 			if (clust_count.find(assignments[i - 1] - 1) == clust_count.end())
	// 			{
	// 				clust_count[assignments[i - 1] - 1] = 0;
	// 			}
	// 			else {
	// 				clust_count[assignments[i - 1] - 1]++;
	// 			}
	// 			//std::cout << "]\n";
	// 		}
	// 		//std::cout << "}\n";

	// 		std::vector<std::vector<int>> super_cluster(nb_super_clusters);
	// 		for (int i = 0; i < super_assignments.size(); i++)
	// 		{
	// 			super_cluster[super_assignments[i]].push_back(i);
	// 		}

	// 		for(int i = 0; i < super_cluster.size(); i++)
	// 		{
	// 			int count = 0;
	// 			for ( int j = 0; j < super_cluster[i].size(); j++)
	// 			{
	// 				count += clust_count[super_cluster[i][j]];
	// 			}
	// 			std::cerr << count << " " << super_cluster[i].size() << "\n";
	// 		}

	// 		Solution best_of_all(&data);
	// 		std::vector<Data> temp_data(cluster.size(), data);

	// 		GlobalSolution global_solution;

	// 		for (int clust_idx = 0; clust_idx < cluster.size(); clust_idx++)
	// 		{
	// 			temp_data[clust_idx].client = {data.client[0]};
	// 			temp_data[clust_idx].nb_clients = 0;

	// 			for (int i = 0; i < cluster[clust_idx].size(); i++)
	// 			{
	// 				temp_data[clust_idx].client.push_back(data.client[cluster[clust_idx][i]]);
	// 				temp_data[clust_idx].client.back().index = cluster[clust_idx][i];
	// 				temp_data[clust_idx].nb_clients++;
	// 				//				std::cerr << clust[i] << " ";
	// 			}
	// 			std::cout << "\n";

	// 			temp_data[clust_idx].ComputeDistanceMatrix();
	// 			temp_data[clust_idx].ComputeRandomCost();

	// 			// Generate initial cluster solution using savings algorithm
	// 			Solution init_clust_sol(&temp_data[clust_idx]);

	// 			for (int i = 1, j = 0; i <= temp_data[clust_idx].nb_clients; i++, j++)
	// 			{
	// 				init_clust_sol.route.push_back({0, i, 0});
	// 				init_clust_sol.route_cost.push_back(temp_data[clust_idx].time_cost[0][i] + temp_data[clust_idx].time_cost[i][0]);
	// 				init_clust_sol.cost += init_clust_sol.route_cost.back();
	// 			}

	// 			RunSavings(init_clust_sol, temp_data[clust_idx]);

	// 			for (int r = 0; r < init_clust_sol.route.size(); r++)
	// 			{
	// 				auto &route = init_clust_sol.route[r];
	// 				best_of_all.route.push_back({0});

	// 				for (int i = 1; i < route.size() - 1; i++)
	// 				{
	// 					best_of_all.route.back().push_back(temp_data[clust_idx].client[route[i]].index);
	// 				}

	// 				best_of_all.route.back().push_back(0);
	// 				best_of_all.cost += init_clust_sol.route_cost[r];
	// 			}

	// 			global_solution.AddPartialSolution(init_clust_sol);

	// 		}

	// 		global_solution.Print();

	// 		//std::cerr << "A\n";
	// 		//best_of_all.ExportToDot("solution0.png");
	// 		//std::cerr << "A\n";
	// 		
	// 		// return 0;
	// 		best_of_all = Solution(&data);

	// 		std::vector<std::vector<std::vector<int>>> clust_route(nb_clusters);
	// 		std::vector<std::vector<double>> clust_route_cost(nb_clusters);

	// 		for (int clust_idx = 0; clust_idx < cluster.size(); clust_idx++)
	// 		{

	// 			std::cerr << clust_idx << "\n";
	// 			// Perform iterated local search on each cluster
	// 			Solution best(&temp_data[clust_idx]);
	// 			best.cost = 999999999;

	// 			SetPool my_pool(&temp_data[clust_idx]);
	// 			IlsRvnd ils_rvnd(&temp_data[clust_idx], &my_pool);

	// 			ils_rvnd.global_solution = &global_solution;
	// 			ils_rvnd.clust_idx = &clust_idx;

	// 			//			auto start = std::chrono::system_clock::now();
	// 			ils_rvnd.Run(1, temp_data[clust_idx].nb_clients + 5 * temp_data[clust_idx].nb_vehicles, temp_data[clust_idx].tolerance);
	// 			// ils_rvnd.Run(temp_data[clust_idx].max_iter_a, temp_data[clust_idx].nb_clients + 5 * temp_data[clust_idx].nb_vehicles, temp_data[clust_idx].tolerance);
	// 			//			auto end = std::chrono::system_clock::now();

	// 			best = ils_rvnd.best;
	// 			ils_rvnd.best.cost = 0;
	// 			//			std::chrono::duration<double> elapsedSeconds = end - start;

	// 			//std::cout << "**********************************************************************************************" << std::endl;
	// 			//std::cout << ils_rvnd.best << std::endl;
	// 			//std::cout << "**********************************************************************************************" << std::endl;
	// 			//std::cout << "AlgTime " << (elapsedSeconds.count()) << " \n\n";

	// 			//// Overwrite BKS in case it was improved
	// 			//if (arg.params.find("-sol") != arg.params.end())
	// 			//{
	// 			//	Solution BKS(&temp_data[clust_idx]);
	// 			//	BKS.ReadSolFile(arg.params["-sol"]);
	// 			//	if (best.cost < BKS.cost)
	// 			//	{
	// 			//		best.WriteSolFile(arg.params["-sol"]);
	// 			//	}
	// 			//}

	// 			clust_route[clust_idx] = std::vector<std::vector<int>> (best.route.size());
	// 			clust_route_cost[clust_idx] = std::vector<double> (best.route.size());

	// 			for (int r = 0; r < best.route.size(); r++)
	// 			{
	// 				// test_cost += best.route_cost[r]; 
	// 				clust_route_cost[clust_idx][r] = best.route_cost[r];

	// 				auto &route = best.route[r];
	// 				best_of_all.route.push_back({0});

	// 				clust_route[clust_idx][r].push_back(0);

	// 				for (int i = 1; i < route.size() - 1; i++)
	// 				{
	// 					best_of_all.route.back().push_back(cluster[clust_idx][route[i] - 1]);
	// 					clust_route[clust_idx][r].push_back(route[i]);
	// 				}

	// 				best_of_all.route.back().push_back(0);
	// 				clust_route[clust_idx][r].push_back(0);
	// 				best_of_all.cost += best.route_cost[r];
	// 			}

	// 			global_solution.partial_solution[clust_idx] = best;
	// 			global_solution.Print();

	// 		}

	// 		best_of_all.ExportToDot("solution1.png");

	// 		// building new global solution by
	// 		// merging clusters into super clusters
	// 		std::vector<Solution> temp_sol_vec = global_solution.partial_solution;
	// 		global_solution.partial_solution.clear();
	// 		temp_data = std::vector<Data>(super_cluster.size(), data);

	// 		std::vector<std::vector<int>> cli_idx_offset(super_cluster.size());
	// 		for (int super_cluster_idx = 0; super_cluster_idx < super_cluster.size(); super_cluster_idx++)
	// 		{
	// 			if (super_cluster[super_cluster_idx].size() < 1)
	// 				continue;

	// 			auto &curr_sup = super_cluster[super_cluster_idx];

	// 			temp_data.push_back(data);
	// 			temp_data[super_cluster_idx].client = {data.client[0]};
	// 			temp_data[super_cluster_idx].nb_clients = 0;
	// 			cli_idx_offset[super_cluster_idx] = {0};

	// 			for (int i = 0; i < curr_sup.size(); i++)
	// 			{
	// 				int offset = cli_idx_offset[super_cluster_idx].back();
	// 				for (int j = 0; j < cluster[curr_sup[i]].size(); j++)
	// 				{
	// 					temp_data[super_cluster_idx].client.push_back(data.client[cluster[curr_sup[i]][j]]);
	// 					temp_data[super_cluster_idx].client.back().index = cluster[curr_sup[i]][j];
	// 					temp_data[super_cluster_idx].nb_clients++;
	// 					offset++;
	// 				}
	// 				cli_idx_offset[super_cluster_idx].push_back(offset);
	// 			}

	// 			temp_data[super_cluster_idx].ComputeDistanceMatrix();
	// 			temp_data[super_cluster_idx].ComputeRandomCost();

	// 			// building initial solution
	// 			// by merging the routes of multiple clusters 

	// 			Solution init_clust_sol(&temp_data[super_cluster_idx]);

	// 			for (int i = 0; i < curr_sup.size(); i++)
	// 			{
	// 				for (int r = 0; r < clust_route[curr_sup[i]].size(); r++)
	// 				{
	// 					init_clust_sol.route.push_back({0});
	// 					init_clust_sol.route_cost.push_back(clust_route_cost[curr_sup[i]][r]);
	// 					init_clust_sol.cost += init_clust_sol.route_cost.back();

	// 					for (int j = 1; j < clust_route[curr_sup[i]][r].size() - 1; j++)
	// 					{
	// 						init_clust_sol.route.back().push_back(cli_idx_offset[super_cluster_idx][i] + clust_route[curr_sup[i]][r][j]);
	// 					}
	// 					init_clust_sol.route.back().push_back(0);

	// 					init_clust_sol.route_load.push_back(0);
	// 					for (int j = 1; j < (int) init_clust_sol.route.back().size() - 1; j++)
	// 					{
	// 						int cli = init_clust_sol.route.back()[j];
	// 						std::cout << init_clust_sol.route.back().size() << " " << cli << "\n";
	// 						init_clust_sol.route_load.back() += temp_data[super_cluster_idx].client[cli].demand;
	// 					}
	// 				}
	// 			}

	// 			init_clust_sol.InitializeNeighbStatus();

	// 			global_solution.AddPartialSolution(init_clust_sol);
	// 		}

	// 		// std::cerr << test_cost << "\n";

	// 		// perform ILS between close clusters
	// 		// (i.e., super-cluster-wise ILS)

	// 		Solution best_of_all_2(&data);

	// 		std::cerr << "Inter cluster phase\n";
	// 		for (int super_cluster_idx = 0; super_cluster_idx < super_cluster.size(); super_cluster_idx++)
	// 		{
	// 			if (super_cluster[super_cluster_idx].size() < 1)
	// 				continue;

	// 			auto &curr_sup = super_cluster[super_cluster_idx];

	// 			//temp_data.push_back(data);
	// 			//temp_data[super_cluster_idx].client = {data.client[0]};
	// 			//temp_data[super_cluster_idx].nb_clients = 0;
	// 			//cli_idx_offset[super_cluster_idx] = {0};

	// 			//for (int i = 0; i < curr_sup.size(); i++)
	// 			//{
	// 			//	int offset = cli_idx_offset[super_cluster_idx].back();
	// 			//	for (int j = 0; j < cluster[curr_sup[i]].size(); j++)
	// 			//	{
	// 			//		temp_data[super_cluster_idx].client.push_back(data.client[cluster[curr_sup[i]][j]]);
	// 			//		temp_data[super_cluster_idx].client.back().index = cluster[curr_sup[i]][j];
	// 			//		temp_data[super_cluster_idx].nb_clients++;
	// 			//		//offset++;
	// 			//	}
	// 			//	cli_idx_offset[super_cluster_idx].push_back(offset);
	// 			//}

	// 			//temp_data[super_cluster_idx].ComputeDistanceMatrix();
	// 			//temp_data[super_cluster_idx].ComputeRandomCost();

	// 			//Solution init_clust_sol(&temp_data[super_cluster_idx]);
	// 			//double test1 = 0;
	// 			//init_clust_sol.route_load = {};

	// 			//for (int i = 0; i < curr_sup.size(); i++)
	// 			//{
	// 			//	for (int r = 0; r < clust_route[curr_sup[i]].size(); r++)
	// 			//	{
	// 			//		init_clust_sol.route.push_back({0});
	// 			//		init_clust_sol.route_cost.push_back(clust_route_cost[curr_sup[i]][r]);
	// 			//		test1 += init_clust_sol.route_cost.back();
	// 			//		init_clust_sol.cost += init_clust_sol.route_cost.back();

	// 			//		for (int j = 1; j < clust_route[curr_sup[i]][r].size() - 1; j++)
	// 			//		{
	// 			//			init_clust_sol.route.back().push_back(cli_idx_offset[super_cluster_idx][i] + clust_route[curr_sup[i]][r][j]);
	// 			//		}
	// 			//		init_clust_sol.route.back().push_back(0);

	// 			//		init_clust_sol.route_load.push_back(0);
	// 			//		for (int j = 1; j < (int) init_clust_sol.route.back().size() - 1; j++)
	// 			//		{
	// 			//			int cli = init_clust_sol.route.back()[j];
	// 			//			std::cout << init_clust_sol.route.back().size() << " " << cli << "\n";
	// 			//			init_clust_sol.route_load.back() += temp_data[super_cluster_idx].client[cli].demand;
	// 			//		}
	// 			//	}
	// 			//}

	// 			// init_clust_sol.InitializeNeighbStatus();

	// 			//double test2 = 0;
	// 			//for (int i = 0; i < init_clust_sol.route.size(); i++)
	// 			//{
	// 			//	for (int j = 0; j < init_clust_sol.route[i].size() - 1; j++)
	// 			//	{
	// 			//		int c1 = init_clust_sol.route[i][j];
	// 			//		int c2 = init_clust_sol.route[i][j + 1];
	// 			//		test2 += temp_data[super_cluster_idx].time_cost[c1][c2];
	// 			//	}
	// 			//}
	// 			//std::cerr << "TEST " << test1 << " " << test2 << '\n';
	// 			//std::cerr << init_clust_sol << "\n";

	// 			//std::cerr << init_clust_sol << "\n";

	// 			Solution init_clust_sol = global_solution.partial_solution[super_cluster_idx];
	// 			init_clust_sol.ExportToDot("solution_temp" + std::to_string(super_cluster_idx) +  ".png");

	// 			IlsRvnd ils_rvnd(&temp_data[super_cluster_idx], NULL, &init_clust_sol);
	// 			ils_rvnd.Run(1, 1000, temp_data[super_cluster_idx].tolerance);

	// 			ils_rvnd.best.ExportToDot("solution_temp_ils" + std::to_string(super_cluster_idx) +  ".png");

	// 			Solution best = ils_rvnd.best;

	// 			for (int r = 0; r < best.route.size(); r++)
	// 			{
	// 				auto &route = best.route[r];
	// 				best_of_all_2.route.push_back({0});

	// 				for (int i = 1; i < route.size() - 1; i++)
	// 				{
	// 					best_of_all_2.route.back().push_back(temp_data[super_cluster_idx].client[best.route[r][i]].index);
	// 				}

	// 				best_of_all_2.route.back().push_back(0);
	// 				best_of_all_2.cost += best.route_cost[r];
	// 			}
	// 		}

	// 		best_of_all_2.ExportToDot("solution2.png");

	// 		// Clarke-wright
	// 		// compute savings list

	// 		RunSavings(best_of_all_2, data);

	// 		best_of_all_2.ExportToDot("solution3.png");

	// 		auto end = std::chrono::system_clock::now();
	// 		std::chrono::duration<double> elapsedSeconds = end - start;

	// 		std::cout << "**********************************************************************************************" << std::endl;
	// 		std::cout << best_of_all_2 << "\n";
	// 		std::cout << "**********************************************************************************************" << std::endl;
	// 		std::cout << "AlgTime " << (elapsedSeconds.count()) << " \n\n";

	// 		for (int i = 0; i < nb_clusters; i++)
	// 		{
	// 			// std::cerr << clust_count[i] << "\n";
	// 		}
	// 	}

	// }
	// else 
	{

		data.ComputeDistanceMatrix();
		data.ComputeRandomCost();
		data.ComputeAuxVectors();

		Solution best_of_all(&data);
		best_of_all.cost = 999999999;

		SetPool my_pool(&data);
		IlsRvnd ils_rvnd(&data, &my_pool);

		auto start = std::chrono::system_clock::now();
		ils_rvnd.Run(data.max_iter_a, data.nb_clients + 5 * data.nb_vehicles, data.tolerance);
		auto end = std::chrono::system_clock::now();

		std::chrono::duration<double> elapsedSeconds = end - start;

		std::cout << "**********************************************************************************************" << std::endl;
		std::cout << ils_rvnd.best << std::endl;
		std::cout << "**********************************************************************************************" << std::endl;
		std::cout << "AlgTime " << (elapsedSeconds.count()) << " \n\n";

		// Overwrite BKS in case it was improved
		if (arg.params.find("-sol") != arg.params.end())
		{
			Solution BKS(&data);
			BKS.ReadSolFile(arg.params["-sol"]);
			if (best_of_all.cost < BKS.cost)
			{
				best_of_all.WriteSolFile(arg.params["-sol"]);
			}
		}

	}


	// best_of_all.ExportToDot("a");
}

bool RunSavings(Solution &best_of_all, Data &data)
{
	std::vector<std::tuple<double, int, int>> saving_list;
	std::vector<int> node_route(data.nb_clients + 1, 0);
	std::vector<int> client_memo(data.nb_clients + 1, 0);
	std::vector<double> load(best_of_all.route.size());

	for (int r1 = 0; r1 < best_of_all.route.size(); r1++)
	{
		if (best_of_all.route[r1].size() <= 2)
			continue;

		load[r1] = 0;
		for (int i = 1; i < best_of_all.route[r1].size() - 1; i++)
		{
			load[r1] += data.client[best_of_all.route[r1][i]].demand;
		}

		int i = best_of_all.route[r1].end()[-2];
		node_route[i] = r1;
		node_route[best_of_all.route[r1][1]] = r1;

		for (int r2 = 0; r2 < best_of_all.route.size(); r2++)
		{
			if (r1 == r2 || best_of_all.route[r2].size() <= 2)
				continue;

			int j = best_of_all.route[r2][1];
			double saving = - data.GetDist(i, 0) - data.GetDist(0, j) + data.GetDist(i, j);

			saving_list.push_back({saving, i, j});
		}
	}

	std::sort(saving_list.begin(), saving_list.end());
	double reduced_cost = 0;

	bool improved = false;

	for (int rr = 0; rr < saving_list.size(); rr++)
	{
		int saving = std::get<0>(saving_list[rr]);

		if (saving >= 0)
			break;

		int i = std::get<1>(saving_list[rr]);
		int j = std::get<2>(saving_list[rr]);
		int r1 = node_route[i];
		int r2 = node_route[j];

		if (r1 == r2 
				|| best_of_all.route[r1].size() <= 2
				|| best_of_all.route[r2].size() <= 2
				|| best_of_all.route[r1].end()[-2] != i
				|| best_of_all.route[r2][1] != j
				|| load[r1] + load[r2] > data.capacity)
		{
			continue;
		}

		//		std::cerr << r1 << " " << r2 << " " << best_of_all.route.size() << "\n";
		best_of_all.route.push_back(std::vector<int>(best_of_all.route[r1].begin(), best_of_all.route[r1].end() - 1));
		best_of_all.route.back().insert(best_of_all.route.back().end(), best_of_all.route[r2].begin() + 1, best_of_all.route[r2].end());

		load.push_back(load[r1] + load[r2]);

		node_route[i] = best_of_all.route.size() - 1;
		node_route[j] = best_of_all.route.size() - 1;
		node_route[best_of_all.route[r1][1]] = best_of_all.route.size() - 1;
		node_route[best_of_all.route[r2].end()[-2]] = best_of_all.route.size() - 1;

		best_of_all.route[r1] = {};
		best_of_all.route[r2] = {};

		reduced_cost += saving;
		best_of_all.cost += saving;
		improved = true;
	}

	for (int r = 0; r < best_of_all.route.size();)
	{
		if (best_of_all.route[r].size() <= 2)
			best_of_all.route.erase(best_of_all.route.begin() + r);
		else ++r;
	}	

	for (auto &route : best_of_all.route)
	{
		best_of_all.route_cost.push_back({0});
		for (int i = 0; i < route.size() - 1; i++)
		{
			best_of_all.route_cost.back() += data.GetDist(route[i], route[i+1]);
		}
	}

	return improved;
}

#else

bool RunSavings(Solution &best_of_all, Data &data);

int main(int argc, char *argv[]) 
{	
	ArgParser arg(argc, argv);


	Data *data = new Data(arg.params);

	data->ComputeDistanceMatrix();
	data->ComputeRandomCost();

	Solution *s = new Solution;
	// *s = GenerateInitialSolution(data, 0, 50);
	s->data = data;

	s->route = 
	{
		{0, 203, 484, 786, 149, 855, 493, 288, 111, 677, 419, 626, 285, 311, 790, 526, 567, 180, 747, 933, 340, 232, 811, 531, 284, 0},
		{0, 141, 186, 969, 724, 561, 297, 22, 728, 745, 219, 960, 881, 474, 175, 101, 294, 94, 406, 41, 588, 972, 501, 80, 660, 301, 0},
		{0, 334, 11, 678, 349, 341, 473, 827, 621, 699, 205, 946, 476, 117, 987, 579, 705, 223, 70, 40, 858, 382, 902, 770, 342, 373, 363, 0},
		{0, 343, 198, 92, 28, 240, 839, 825, 322, 27, 837, 197, 792, 808, 346, 402, 802, 104, 298, 48, 245, 293, 491, 0},
		{0, 637, 68, 411, 953, 554, 982, 445, 932, 880, 7, 176, 399, 351, 594, 510, 639, 534, 605, 688, 357, 409, 874, 945, 207, 0},
		{0, 313, 907, 681, 685, 608, 79, 721, 598, 746, 920, 602, 569, 251, 192, 943, 488, 744, 864, 513, 318, 983, 99, 830, 620, 0},
		{0, 638, 740, 446, 270, 897, 506, 398, 941, 102, 258, 154, 796, 936, 627, 725, 37, 183, 771, 134, 229, 652, 444, 145, 0},
		{0, 690, 307, 788, 635, 148, 614, 360, 872, 347, 516, 465, 644, 686, 572, 661, 846, 215, 190, 697, 773, 429, 737, 0},
		{0, 39, 379, 116, 966, 14, 459, 671, 974, 400, 633, 582, 227, 581, 859, 785, 443, 731, 659, 43, 892, 378, 188, 82, 0},
		{0, 159, 575, 308, 924, 87, 591, 405, 995, 576, 20, 388, 869, 777, 286, 544, 83, 393, 795, 173, 551, 893, 136, 991, 910, 0},
		{0, 457, 179, 372, 955, 146, 734, 934, 489, 139, 957, 206, 431, 774, 959, 242, 290, 467, 938, 417, 51, 613, 86, 913, 819, 0},
		{0, 704, 249, 212, 333, 482, 185, 512, 401, 380, 618, 326, 319, 631, 271, 666, 878, 849, 530, 389, 556, 646, 78, 989, 0},
		{0, 44, 345, 868, 693, 165, 89, 509, 279, 133, 980, 452, 283, 812, 201, 508, 871, 45, 306, 248, 266, 604, 356, 386, 0},
		{0, 218, 309, 573, 782, 269, 155, 100, 885, 152, 667, 169, 940, 178, 36, 470, 997, 224, 715, 69, 636, 750, 807, 824, 0},
		{0, 277, 52, 262, 752, 861, 898, 392, 246, 814, 779, 726, 566, 428, 425, 273, 727, 171, 261, 663, 708, 596, 103, 984, 578, 884, 0},
		{0, 540, 538, 96, 815, 462, 370, 763, 793, 5, 606, 122, 164, 695, 63, 914, 487, 150, 870, 826, 789, 736, 230, 0},
		{0, 656, 755, 330, 759, 135, 160, 904, 259, 263, 821, 475, 838, 138, 574, 939, 72, 310, 928, 733, 831, 545, 414, 954, 0},
		{0, 742, 424, 413, 202, 427, 906, 329, 34, 776, 714, 2, 979, 649, 909, 253, 365, 217, 657, 970, 889, 439, 391, 49, 590, 0},
		{0, 801, 234, 236, 836, 300, 321, 485, 852, 623, 479, 899, 723, 316, 354, 433, 237, 816, 822, 187, 958, 515, 876, 748, 0},
		{0, 517, 471, 95, 296, 6, 895, 654, 3, 161, 973, 756, 901, 555, 599, 216, 520, 629, 91, 32, 50, 454, 645, 761, 670, 0},
		{0, 565, 214, 937, 875, 794, 680, 144, 580, 105, 642, 805, 407, 767, 931, 570, 533, 289, 650, 132, 600, 231, 865, 320, 735, 988, 0},
		{0, 762, 511, 922, 609, 684, 182, 563, 890, 860, 121, 961, 550, 387, 189, 810, 466, 435, 769, 272, 543, 985, 764, 156, 0},
		{0, 903, 766, 921, 524, 844, 804, 702, 942, 944, 355, 634, 9, 361, 15, 611, 900, 935, 151, 696, 595, 120, 369, 505, 0},
		{0, 492, 710, 220, 282, 74, 707, 722, 123, 21, 42, 919, 16, 274, 562, 325, 549, 107, 632, 753, 993, 775, 25, 0},
		{0, 999, 856, 371, 226, 911, 615, 292, 896, 674, 641, 529, 539, 711, 257, 950, 628, 317, 800, 243, 783, 55, 0},
		{0, 976, 851, 191, 35, 823, 208, 423, 963, 19, 420, 967, 239, 235, 162, 477, 835, 546, 410, 651, 926, 883, 324, 0},
		{0, 295, 478, 536, 170, 275, 97, 268, 323, 336, 610, 381, 368, 147, 67, 77, 679, 26, 559, 497, 716, 675, 250, 174, 0},
		{0, 720, 481, 542, 998, 64, 809, 157, 377, 463, 490, 535, 918, 929, 142, 374, 328, 867, 818, 619, 168, 964, 522, 0},
		{0, 857, 854, 1, 244, 359, 447, 847, 541, 181, 362, 780, 689, 287, 421, 385, 344, 828, 975, 158, 304, 834, 278, 0},
		{0, 177, 416, 441, 577, 495, 422, 701, 238, 778, 843, 981, 923, 442, 781, 866, 882, 60, 586, 456, 394, 622, 367, 119, 0},
		{0, 440, 499, 528, 648, 719, 583, 739, 927, 952, 194, 460, 264, 560, 61, 390, 703, 655, 850, 547, 640, 280, 494, 211, 840, 0},
		{0, 585, 709, 765, 129, 630, 364, 222, 523, 994, 228, 558, 47, 820, 33, 968, 437, 568, 947, 163, 597, 196, 468, 455, 172, 464, 0},
		{0, 537, 557, 643, 143, 167, 58, 88, 85, 62, 93, 817, 694, 784, 905, 607, 571, 915, 806, 110, 760, 453, 469, 0},
		{0, 592, 153, 255, 751, 130, 65, 730, 75, 432, 877, 888, 917, 353, 949, 53, 81, 842, 384, 204, 140, 59, 754, 0},
		{0, 925, 616, 717, 480, 738, 396, 990, 305, 128, 327, 54, 691, 886, 76, 302, 500, 438, 833, 548, 647, 706, 114, 213, 221, 0},
		{0, 461, 199, 315, 948, 109, 887, 519, 841, 57, 73, 449, 338, 848, 584, 193, 682, 131, 339, 741, 458, 486, 553, 124, 0},
		{0, 502, 265, 397, 593, 879, 418, 668, 247, 625, 853, 672, 962, 312, 1000, 200, 451, 797, 256, 71, 337, 31, 0},
		{0, 518, 713, 624, 676, 376, 412, 375, 971, 732, 472, 552, 332, 587, 125, 403, 665, 862, 589, 254, 758, 832, 241, 653, 267, 366, 0},
		{0, 299, 712, 395, 564, 525, 38, 603, 965, 718, 692, 209, 799, 450, 118, 514, 291, 430, 434, 90, 23, 108, 617, 0},
		{0, 8, 683, 436, 195, 210, 483, 956, 803, 845, 112, 992, 24, 303, 498, 757, 352, 56, 84, 863, 10, 916, 29, 115, 507, 0},
		{0, 358, 813, 17, 404, 729, 787, 46, 496, 106, 612, 891, 281, 986, 4, 276, 184, 166, 951, 791, 521, 996, 137, 127, 829, 348, 0},
		{0, 331, 448, 912, 415, 503, 687, 66, 798, 772, 113, 743, 314, 30, 426, 252, 700, 504, 527, 225, 233, 749, 669, 977, 383, 126, 0},
		{0, 98, 894, 673, 930, 350, 601, 260, 13, 908, 12, 978, 532, 18, 873, 335, 662, 768, 698, 664, 408, 658, 0}
	};


	//for(int i = 0; i < data->time_cost.size(); i++)
	//{
	//	for (int j = 0; j < data->time_cost[i].size(); j++)
	//	{
	//		std::cout << data->time_cost[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}


	s->route_cost.resize(s->route.size());
	s->route_load.resize(s->route.size());

	s->InitializeNeighbStatus();

	ADS* ads = new ADS(data);
	ads->SetSolution(s);
	ads->SetPenalizedCost();
	std::cout << *s << "\n";

	auto start = std::chrono::system_clock::now();
	LocalSearch* local_search = new LocalSearch(s, ads);
	IlsStats *stats = new IlsStats;
	local_search->stats = stats;
	local_search->Run();
	// local_search->SearchRelocate1Complete();

	ads->SetPenalizedCost();
	std::cout << *s << "\n";
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsedSeconds = end - start;
	std::cout << elapsedSeconds.count() << "\n";

	delete s;
	delete ads;
	delete data;
	delete local_search;
	delete stats;
	//	delete perturbation;
} 
#endif
