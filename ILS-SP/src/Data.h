#ifndef DATA_H
#define DATA_H

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
#include <map>

#define EPSILON 1e-5

static int best_nb_vehi;

struct Vertex
{
	int index;
	double x = 0;
	double y = 0;
	double demand = 0;
	double demand_deviation = 0;
	double serv_time = 0;
	double ready_time = 0;
	double due_time = 0;
};

struct Data
{
	/** Algorithm parameters */
	int seed;
	int granular_size = 20;
	double penalty_load = 10000;
	double penalty_tw = 10000;
	double t_dev_a = 0.05;
	double tolerance = 0.005;
	int max_iter_a = 50;
	int max_iter_ils_a = 1000;
	double t_dev_b = 0.005;
	int max_iter_b = 25;
	int max_iter_ils_b = 200;
	double max_iter_ils;
	double maxRootGap = 0.02;
	double max_sp_time = 60;
	int A = 11;
	int N = 150;
	bool use_acceptance_criterion = false;
	bool use_swap_star = false;

	/** Output info */
	std::string dot_file_path; // .dot file for graphviz Plotting
	std::string sol_file_path; // .sol file in standard CVRPlib format

	/** Problem data */
	std::string name;
	std::string variant;
	int dist_type = 1;
	int nb_clients;
	int lb_vehicles;
	int nb_vehicles;
	double capacity;
	bool limit_nb_vehicles = false;
	std::vector<Vertex> client;				   // set of clients + depot
	std::vector<std::vector<double>> time_cost; // travel cost
	std::vector<std::vector<double>> time_deviation;
	int budget_time = 1;
	int budget_demand = 5;
	double alpha_time = 0.1;
	double alpha_demand = 0.0;
	int min_cap, max_cap;
	double dem_dev_factor = 0.3;
	double avg_route_len_factor = 0.75;
	double cap_weight = 0.3;
	std::vector<double> demand;
	std::vector<double> demand_deviation;

	/** Preprocessed data */
	std::vector<std::vector<int>> correlated_vertices;
	std::vector<std::vector<double>> random_cost;

	// ins_before[i][j] = 1 if client i may come before 
	// client j in a feasible solution, 0 otherwise.
	std::vector<std::vector<int>> ins_before; 

	/** Instance reading/generating methods */
	Data(std::map<std::string, std::string> &params);
	Data(int argc, char **argv);
	Data(std::string instancePath, int variant, int seed);

	void LoadInstanceCVRP(std::string &fileName);
	void LoadInstancePRP(int argc, char **argv);
	void LoadInstanceVRPTW(std::string &fileName);
	void ComputeDistanceMatrix();
	void ComputeAuxVectors();
	int GetDist(int i, int j) { return std::floor(std::sqrt(std::pow(client[i].x - client[j].x, 2) + std::pow(client[i].y - client[j].y, 2)) + 0.5); }
	void CalculateProximity();
	void ComputeRVRPTWParameters();
	void ComputeRandomCost();
	void PreProcessFsr();

	int FindMinFeasibleCapacity(int nb_vehicles); // used specifically for CVRP instances by Pessoa et al.
	int CalcMinNbVehicles();
	void ComputeRCVRPParameters();

	void GetName(std::string &instanceName);
};

#endif
