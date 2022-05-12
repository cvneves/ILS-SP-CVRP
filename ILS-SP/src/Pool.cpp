#include "Pool.h"

double SetPool::GetHash(std::vector<int> &sequence)
{
	double hash = 0;
	for (int i = 0; i < sequence.size() - 1; i++)
	{
		hash += data->random_cost[sequence[i]][sequence[i + 1]];
	}
	// std::cout << hash  << "\n";
	return hash;
}

void SetPool::Update(Solution *s, bool is_permanent)
{

	for(int r = 0; r < s->route.size(); r++)
	{
		// if (s->route[r].size() > 2 && !is_permanent && s->route_stats[r].solution_gap < 1.5)
		if (s->route[r].size() > 2 && !is_permanent)
		{
			bool found = false;
			for (int i = 0; i < temp.size(); i++)
			{
				if (s->route[r][1] == temp[i].sequence[1] // 
						&& s->route[r].end()[-2] == temp[i].sequence.end()[-2]  //
						&& s->route[r][s->route[r].size() / 2] == temp[i].sequence[temp[i].sequence.size() / 2] //
						&& s->route_load[r] == temp[i].load //
						&& s->route_cost[r] == temp[i].cost //
						&& s->route[r].size() == temp[i].sequence.size()) //
				{
					{
						found = true;
						break;
					}
				}
			}

			if (!found)
			{
					temp.push_back({s->route[r], s->route_cost[r], s->route_load[r], s->route_stats[r]});
			}
		}

		else if (s->route[r].size() > 2){
		{
			bool found = false;
			for (int i = 0; i < permanent.size(); i++)
			{
				if (s->route[r][1] == permanent[i].sequence[1] // 
						&& s->route[r].end()[-2] == permanent[i].sequence.end()[-2]  //
						&& s->route[r][s->route[r].size() / 2] == permanent[i].sequence[permanent[i].sequence.size() / 2] //
						&& s->route_load[r] == permanent[i].load //
						&& s->route_cost[r] == permanent[i].cost //
						&& s->route[r].size() == permanent[i].sequence.size()) //
				{
					{
						found = true;
						break;
					}
				}
			}

			if (!found)
			{
					permanent.push_back({s->route[r], s->route_cost[r], s->route_load[r], s->route_stats[r]});
			}
		}

					//permanent.push_back({s->route[r], s->route_cost[r], s->route_load[r]});
		}
	}

	

	// for (int r = 0; r < s->route.size(); r++)
	// {
	// 	if (s->route[r].size() <= 2 || s->route_stats[r].solution_gap > gap_tol)
	// 		continue;

	// 	double hash = GetHash(s->route[r]);

	// 	bool found_in_perm = p_key.find(hash) != p_key.end();

	// 	if (is_permanent && !found_in_perm)
	// 	{
	// 		permanent.push_back({s->route[r], s->route_cost[r], s->route_load[r], s->route_stats[r]});
	// 		p_key.insert(hash);
	// 	}

	// 	if (!is_permanent && !found_in_perm && (np_key.find(hash) == np_key.end()))
	// 	{
	// 		temp.push_back({s->route[r], s->route_cost[r], s->route_load[r], s->route_stats[r]});
	// 		np_key.insert(hash);
	// 	}
	// }
}

void SetPool::ForceTemp(Solution *s)
{
	for (int r = 0; r < s->route.size(); r++)
	{
		if (s->route[r].size() <= 2)
			continue;

		double hash = GetHash(s->route[r]);
		np_key.insert(hash);
		temp.push_back({s->route[r], s->route_cost[r], s->route_load[r], s->route_stats[r]});
		// temp.push_back({s->route[r], s->route_cost[r], s->route_load[r], solution_gap, nb_it_ils, nb_restart});
	}
}

void SetPool::RemoveTempRoutes() { temp.clear(); np_key.clear();}
void SetPool::GapClear(double best_sol_cost)
{
	int nb_removed = 0;

	std::cout << "Cleaning pool by gap..." << std::endl;

	for (auto it = temp.begin(); it != temp.end();)
	{
		if (((double) 100.0 * (it->stats.solution_cost - best_sol_cost) / best_sol_cost) > gap_tol)
		{
			double hash = GetHash(it->sequence);
			np_key.erase(hash);

			nb_removed++;
			it = temp.erase(it);
		}
		else
		{
			++it;
		}
	}

	std::cout << nb_removed << " routes removed from pool" << std::endl;
}

RouteInfo *SetPool::GetPermRoute(int p) { return &permanent[p]; }
RouteInfo *SetPool::GetTempRoute(int p) { return &temp[p]; }
int SetPool::GetNbPermPath() { return permanent.size(); }
int SetPool::GetNbNonPermPath() { return temp.size(); }
