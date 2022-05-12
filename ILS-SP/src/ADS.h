#ifndef ADS_H
#define ADS_H

#include <vector>
#include "MoveMemo.h"
#include "RouteADS.h"
#include "Solution.h"

enum SMD_STATUS
{
	UNVISITED = -1, UNIMPROVED = 0,	IMPROVED = 1
};

struct SMD
{
	int r1; 
	int r2;
	int r1_pos1;
	int r1_pos2;
	int r2_pos1;
	int r2_pos2;
	int move_type;

	SMD_STATUS status;

	double cost;

	SMD Reverse()
	{
		int mt = move_type;
		if (mt==3) mt=2;
		else if (mt==2) mt=3;
		return {
			r2,
			r1,
			r2_pos1,
			r2_pos2,
			r1_pos1,
			r1_pos2,
			mt,
			status,
			cost
		};
	}
};

struct ADS
{
	Data *data;
	std::vector<RouteADS*> routeADS;
	Solution *curr_solution;
	int used_vehicles;
	int contador1 = 0;
	int contador2 = 0;

	// Solution memoization
	std::unordered_map<double, int> solution_memo;
	double GetCurrSolutionHash();
	bool FindAndPush(int iter); // adds Solution to the list and returns true if it was already there

	// Movement memoization
	std::map<std::pair<double, double>, std::vector<SMD>> inter_mov_memo;
	std::map<double, std::vector<SMD>> intra_mov_memo;
	std::pair<double, double> GetRoutePairHash(int r1, int r2);
	SMD* getSMD(int r1, int r2, int neighborhood);
	SMD* getSMD(int r, int neighborhood);
	void push(int r1, int r2, int neighborhood);
	void push(int r, int neighborhood);


	void SetSolution(Solution *s);
	void SetPenalizedCost();
	ADS(Data *data) : data(data) {}
	~ADS();
};

#endif
