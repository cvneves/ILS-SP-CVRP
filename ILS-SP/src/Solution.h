#ifndef SOLUTION_H
#define SOLUTION_H

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
#include <sstream>
#include <set>

#include "Data.h"
#include "RouteStats.h"

struct Solution
{
	Data *data;
	std::vector<std::vector<int>> route;
	std::vector<double> route_cost;
	std::vector<double> route_load;

	double cost = std::numeric_limits<double>::infinity();
	int nb_vehicles;

	std::vector<RouteStats> route_stats; 

	// ADS *ads;

	Solution() {}
	Solution(Data *data) { this->data = data; }

	// Neighborhood status
	std::vector<std::vector<std::vector<int>>> inter_neighb_status;
	std::vector<std::vector<int>> intra_neighb_status;

	std::set<std::pair<int, int>> GetArcSet();
	int GetDistance(Solution *s);

	void InitializeNeighbStatus();
	void SetInterNeighbStatusAllTrue(int r1, int r2);
	void SetInterNeighbStatusAllTrue(int r1);
	void SetInterNeighbStatus(int r1, int r2, int neighb, bool status);
	bool GetInterNeighbStatus(int r1, int r2, int neighb);
	void SetIntraNeighbStatus(int r, int neighb, bool status) {intra_neighb_status[r][neighb] = status;}
	bool GetIntraNeighbStatus(int r, int neighb) {return intra_neighb_status[r][neighb];}
	void SetIntraNeighbStatusAllTrue(int r);

	void PrintPartial(int &route_counter);
	void ParallelInsertion();
	void SequentialInsertion();
	void ComputeCost();
	void ExportToDot(std::string dot_file_path);
	void ReadSolFile(std::string sol_file_path);
	void WriteSolFile(std::string sol_file_path);
};

std::ostream &operator<<(std::ostream &os, Solution const &s);

struct GlobalSolution
{
	std::vector<Solution> partial_solution;

	void Print() 
	{
		int route_count = 1;
		for (auto &sol : partial_solution)
		{
			sol.PrintPartial(route_count);
		}
		std::cout << "Cost 0\n";
	}

	void AddPartialSolution(Solution &solution)
	{
		partial_solution.push_back(solution);
	}

	void DeletePartialSolutions(std::vector<int> &index_vec)
	{
		std::sort(index_vec.begin(), index_vec.end());
		while(!index_vec.empty())
		{
			int index = index_vec.back();
			index_vec.pop_back();
			partial_solution.erase(partial_solution.begin() + index);
		}
	}
};

#endif
