#ifndef PERTURBATION_H
#define PERTURBATION_H

#include "Solution.h"
#include "ADS.h"

struct Perturbation
{
	Data *data;
	Solution *curr_solution;
	ADS *ads;
	bool perturbed;
	int last_perturb = -1;

	void MultipleShift();
	void MultipleSwap11();
	void DoubleBridge(int r);

	public:
	std::vector<int> removed_clients;
	std::vector<std::pair<int, int>> adjacent;
	void ComputeAdjacentNodes();
	void RemoveClient(int route, int pos);
	void InsertClient(int route, int pos, int client);
	void ConcentricRemoval(int num_clients);
	void RandomRemoval(int num_clients);
	void SequenceRemoval(int num_clients);
	void InsertionByDistance(int curr_client);

	Perturbation(Solution *curr_solution, ADS *ads);
	void Run(int degree);
};

#endif
