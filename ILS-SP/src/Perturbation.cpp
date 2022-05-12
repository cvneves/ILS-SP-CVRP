#include "Perturbation.h"

Perturbation::Perturbation(Solution* curr_solution, ADS* ads) : data(curr_solution->data), curr_solution(curr_solution), ads(ads)
{

}

void Perturbation::Run(int degree)
{
	perturbed = false;
	if (ads->used_vehicles == 1)
	{
		for (int r = 0; r < curr_solution->route.size(); r++)
		{
			if (curr_solution->route[r].size() > 2)
			{
				DoubleBridge(r);
				break;
			}
		}
		return;
	}

	// RandomRemoval(degree);	
 	// for (auto &i : removed_clients)
 	// 	InsertionByDistance(i);

	// TODO: TIREI UMA PERTURBACAO PARA TESTAR
	// if (rand() % 2)
	while (!perturbed)
	{
		// if (rand() % 2)
		if (rand() % 2)
		{
			MultipleSwap11();
			last_perturb = 0;
		}
		else
		{
			MultipleShift();
			last_perturb = 1;
		}
	}
	// else {
	// 	MultipleShift();
	// }
}

void Perturbation::DoubleBridge(int r)
{
	int c1, c2, c3, c4;

	if (curr_solution->route[r].size() > 8) {
		c1 = 1 + rand() % ((int) curr_solution->route[r].size() - 3);
		c2 = c1 + 1;

		c3 = c1;
		c4 = c1;
		while (c3 == c1 || c3 == c2 || c4 == c1 || c4 == c2) {
			c3 = 1 + rand() % ((int) curr_solution->route[r].size() - 3);
			c4 = c3 + 1;
		}

		std::swap(curr_solution->route[r][c1], curr_solution->route[r][c3]);
		std::swap(curr_solution->route[r][c2], curr_solution->route[r][c4]);

		ads->routeADS[r]->Update();
	}
}

void Perturbation::MultipleSwap11()
{
	std::vector<bool> modified_routes(curr_solution->route.size(), false);
	int k = rand() % 2 + 1;
	int a = 0;
	perturbed = false;

	while (a < k)
	{
		double load_r1, load_r2;
		int r1, r2, k, l;

		r1 = rand() % curr_solution->route.size();
		r2 = r1;
		while (r1 == r2)
			r2 = rand() % curr_solution->route.size();

		if (curr_solution->route[r1].size() > 2 && curr_solution->route[r2].size() > 2)
		{
			load_r1 = load_r2 = 0;
			k = rand() % (curr_solution->route[r1].size() - 2) + 1;
			l = rand() % (curr_solution->route[r2].size() - 2) + 1;

			for (int i = 0; i < curr_solution->route[r1].size(); i++)
				load_r1 += data->client[curr_solution->route[r1][i]].demand;
			for (int i = 0; i < curr_solution->route[r2].size(); i++)
				load_r2 += data->client[curr_solution->route[r2][i]].demand;

			load_r1 += (-data->client[curr_solution->route[r1][k]].demand + data->client[curr_solution->route[r2][l]].demand);
			load_r2 -= (-data->client[curr_solution->route[r1][k]].demand + data->client[curr_solution->route[r2][l]].demand);

			if (load_r1 <= data->capacity && load_r2 <= data->capacity)
			{
				std::swap(curr_solution->route[r1][k], curr_solution->route[r2][l]);

				modified_routes[r1] = true;
				modified_routes[r2] = true;
				perturbed = true;
			}
		}

		a++;
	}

	for (int r = 0; r < modified_routes.size(); r++)
	{
		if (modified_routes[r])
		{
			// std::cout << r << "\n";
			// TODO: (MAYBE) port all Update functions to ADS class
			ads->routeADS[r]->Update();
			curr_solution->SetIntraNeighbStatusAllTrue(r);
			curr_solution->SetInterNeighbStatusAllTrue(r);
		}
	}

	return;

	// ComputePenalizedCost(curr_solution);
}

void Perturbation::MultipleShift()
{
	std::vector<bool> modified_routes(curr_solution->route.size(), false);
	int num_shifts = rand() % 2 + 1;
	int a = 0;
	perturbed = false;

	while (a++ < num_shifts)
	{
		double load_r1 = 0;
		double load_r2 = 0;
		int r1, r2, k, l, k_index, l_index;

		r1 = rand() % curr_solution->route.size();
		r2 = r1;
		while (r2 == r1)
			r2 = rand() % curr_solution->route.size();

		int r1_size = ads->routeADS[r1]->length;
		int r2_size = ads->routeADS[r2]->length;

		if (r1_size > 3 && r2_size > 3)
		{
			k_index = 1 + rand() % (r1_size - 2);
			l_index = 1 + rand() % (r2_size - 2);
			k = curr_solution->route[r1][k_index];
			l = curr_solution->route[r2][l_index];

			load_r1 = 0;
			load_r2 = 0;
			for (int i = 0; i < r1_size; i++)
				load_r1 += data->client[curr_solution->route[r1][i]].demand;
			for (int i = 0; i < r2_size; i++)
				load_r2 += data->client[curr_solution->route[r2][i]].demand;
			load_r1 -= data->client[k].demand - data->client[l].demand;
			load_r2 += data->client[k].demand - data->client[l].demand;

			if (load_r1 + EPSILON <= data->capacity && load_r2 + EPSILON <= data->capacity)
			{
				perturbed = true;

				int r1_index = 1 + rand() % (curr_solution->route[r1].size() - 2);
				int r2_index = 1 + rand() % (curr_solution->route[r2].size() - 2);
				
				curr_solution->route[r1].erase(curr_solution->route[r1].begin() + k_index);
				curr_solution->route[r2].erase(curr_solution->route[r2].begin() + l_index);

				curr_solution->route[r1].insert(curr_solution->route[r1].begin() + r1_index, l);
				curr_solution->route[r2].insert(curr_solution->route[r2].begin() + r2_index, k);

				modified_routes[r1] = true;
				modified_routes[r2] = true;
			}
		}
	}

	for (int r = 0; r < modified_routes.size(); r++)
	{
		if (modified_routes[r])
		{
			ads->routeADS[r]->Update();
			curr_solution->SetIntraNeighbStatusAllTrue(r);
			curr_solution->SetInterNeighbStatusAllTrue(r);
		}
	}

	// ComputePenalizedCost(curr_solution);
}

void Perturbation::ComputeAdjacentNodes()
{
	adjacent.assign(data->nb_clients + 1, std::pair<int,int>());

	for (int r = 0; r < curr_solution->route.size(); r++)
	{
		for (int pos = 1; pos < curr_solution->route[r].size() - 1; pos++)
		{
			int client = curr_solution->route[r][pos];
			int prev = curr_solution->route[r][pos - 1];
			int next = curr_solution->route[r][pos + 1];
			adjacent[client] = {prev, next};
		}
	}
}

void Perturbation::RemoveClient(int r, int pos)
{
	curr_solution->route[r].erase(curr_solution->route[r].begin() + pos);
}

void Perturbation::InsertClient(int r, int pos, int client)
{
	curr_solution->route[r].insert(
		curr_solution->route[r].begin() + pos + 1,
		client
	);
}

void Perturbation::RandomRemoval(int num_clients)
{
	removed_clients.clear();
	ComputeAdjacentNodes();

	for (int k = 0; k < num_clients; k++)
	{
		int r;
		do {
			r = rand() % curr_solution->route.size();
		} while (curr_solution->route[r].size() <= 2);
		int pos = 1 + rand() % (curr_solution->route[r].size() - 2);

		int client = curr_solution->route[r][pos];
		removed_clients.push_back(client);
		RemoveClient(r, pos);	
	}
}

void Perturbation::InsertionByDistance(int curr_client)
{
	double best_cost = std::numeric_limits<double>::infinity();
	int best_pos = 0;
	int best_route = 0;

	for (int r = 0; r < curr_solution->route.size(); r++)
	{
		for (int pos = 1; pos < curr_solution->route[r].size() - 1; pos++)
		{
			int client = curr_solution->route[r][pos];
			int prev = curr_solution->route[r][pos - 1];
			int next = curr_solution->route[r][pos + 1];
			double cost = data->time_cost[curr_solution->route[r][pos]][curr_client];
			if (client != adjacent[curr_client].first && client != adjacent[curr_client].second
					&& next != adjacent[curr_client].first && next != adjacent[curr_client].second
					&& cost < best_cost
				)
			{
				best_cost = cost;	
				best_pos = pos;
				best_route = r;
			}
		}
	}

	InsertClient(best_route, best_pos, curr_client);
}
