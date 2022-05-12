#include "Solution.h"
#include "ADS.h"

std::ostream &operator<<(std::ostream &os, Solution const &s)
{
	int route_count = 0;
	for (int i = 0; i < s.route.size(); i++)
	{
		if (s.route[i].size() <= 2)
		    continue;

		os << "Route #" << route_count++ + 1 << ": ";

		for (int j = 1; j < s.route[i].size() - 1; j++)
		{
			os << s.route[i][j] << " ";
		}
		os << "\n";
	}
	os << "Cost " << s.cost;
	return os;
}

void Solution::PrintPartial(int &route_counter)
{
	for (int i = 0; i < route.size(); i++)
	{
		if (route[i].size() <= 2)
		    continue;

		std::cout << "Route #" << route_counter++ << ": ";

		for (int j = 1; j < route[i].size() - 1; j++)
		{
			std::cout << data->client[route[i][j]].index << " ";
		}
		std::cout << "\n";
	}
}

void Solution::ReadSolFile(std::string sol_file_path)
{
	std::ifstream infile(sol_file_path);

	if (infile.fail())
	{
		cost = 9999999999999;
		return;
	}

	std::string line;
	std::string substr;

	this->route = std::vector<std::vector<int>>(data->nb_vehicles, {0, 0});
	this->route_cost = std::vector<double>(data->nb_vehicles, 0);
	this->route_load = std::vector<double>(data->nb_vehicles, 0);
	int route_count = 0;
	while (std::getline(infile, line))
	{
		std::stringstream linestream(line);
		linestream >> substr;
		std::vector<int> v = {0};

		if (substr == "Route")
		{
			linestream >> substr;
			while (linestream >> substr)
			{
				// std::cout << substr << " ";
				v.push_back(std::atoi(substr.c_str()));
			}
			v.push_back(0);
			this->route[route_count++] = v;

			// std::cout << "\n";
		}

		if (substr == "Cost")
		{
			linestream >> substr;
			this->cost = std::stod(substr.c_str());
		}

		// while (linestream >> substr)
		// {

		//     std::cout << substr << "\n";
		// }
	}

	infile.close();
}

void Solution::WriteSolFile(std::string sol_file_path)
{
	std::ofstream outfile(sol_file_path);
	outfile << *this;
	outfile.close();
}

void Solution::ExportToDot(std::string image_path)
{
	std::ofstream outfile("solution.dot");

	outfile << "digraph G {\n";
	double scale_factor = 500.0;

	outfile << "v" << 0 << " [shape = point,fill=black,pos=\"" << data->client[0].x / scale_factor << "," << data->client[0].y / scale_factor << "!\"]\n";
	for (int i = 1; i <= data->nb_clients; i++)
	{
		outfile << "v" << i << " [shape = point,pos=\"" << data->client[i].x / scale_factor << "," << data->client[i].y / scale_factor << "!\"]\n";
	}

	std::vector<std::string> colors = {"blue",
		"chartreuse",
		"cadetblue",
		"chocolate",
		"darkorchid",
		"deeppink",
		"coral",
		"mediumpurple",
		"red"};

	for (int r = 0; r < route.size(); r++)
	{
		if (route[r].size() == 2)
			continue;
		std::string color = colors[rand() % colors.size()] + std::to_string(1 + rand() % 3);
		for (int i = 0; i < route[r].size() - 1; i++)
		{
			if (i == 0 || i == route[r].size() - 2)
			{
				continue;
				outfile << "v" << route[r][i] << " -> "
					<< "v" << route[r][i + 1] << "[style=dashed,arrowhead=none,penwidth=1,color=" << color << "]"
					<< "\n";
			}
			else
			{
				outfile << "v" << route[r][i] << " -> "
					<< "v" << route[r][i + 1] << "[arrowhead=none,penwidth=5,color=" << color << "]"
					<< "\n";
			}
		}
	}

	outfile << "}";
	outfile.close();

	system(("neato -Tpng solution.dot -o " + image_path).c_str());
}

void Solution::ParallelInsertion()
{
	/** Initializing the candidate list CL with all clients */

	std::vector<int> CL(data->nb_clients);

	for (int i = 1; i <= data->nb_clients; i++)
		CL[i - 1] = i;

	std::random_shuffle(CL.begin(), CL.end());

	route = std::vector<std::vector<int>>(this->nb_vehicles);
	route_cost = std::vector<double>(this->nb_vehicles);
	std::vector<double> load(this->nb_vehicles, 0.0);

	/** Assigning a single client from CL to each route randomly */
	for (int r = 0; r < this->nb_vehicles; r++)
	{
		if (CL.empty())
		{
			route[r] = {0, 0};
			continue;
		}
		route[r] = {0, CL.back(), 0};
		load[r] = data->client[CL.back()].demand;
		CL.pop_back();
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
			for (int r = 0; r < this->nb_vehicles; r++)
			{
				for (int pos = 0; pos < route[r].size() - 1; pos++)
				{
					double insertionCost;

					if (cheapestInsertion)
					{
						insertionCost = data->time_cost[route[r][pos]][*p] + data->time_cost[*p][route[r][pos + 1]] - data->time_cost[route[r][pos]][route[r][pos + 1]] -
							gamma * (data->time_cost[0][*p] + data->time_cost[*p][0]);
					}
					else
					{
						insertionCost = data->time_cost[route[r][pos]][*p];
					}

					if (insertionCost < bestInsertionCost)
					{
						double totalDemand = load[r] + data->client[*p].demand;
						std::vector<double> sorted_deviations(route[r].size() - 1);
						for (int i = 1; i < route[r].size() - 1; i++)
						{
							sorted_deviations[i - 1] = data->client[route[r][i]].demand_deviation;
							// totalDemand += data->client[route[r][i]].demand;
						}
						sorted_deviations.back() = data->client[route[r][pos]].demand_deviation;

						if (data->budget_demand > 0)
						{
							if (route[r].size() - 2 >= data->budget_demand)
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
			route[bestroute].insert(route[bestroute].begin() + bestPosition, *bestNode);
			load[bestroute] += data->client[(*bestNode)].demand;
			CL.erase(bestNode);
		}
		else
		{
			forceFeasible = false;
		}
	}

}

void Solution::SequentialInsertion()
{
	/** Initializing the candidate list CL with all clients */

	std::vector<int> CL(data->nb_clients);

	for (int i = 1; i <= data->nb_clients; i++)
		CL[i - 1] = i;

	std::random_shuffle(CL.begin(), CL.end());

	route = std::vector<std::vector<int>>(this->nb_vehicles);
	route_cost = std::vector<double>(this->nb_vehicles);
	std::vector<double> load(this->nb_vehicles, 0.0);

	/** Assigning a single client from CL to each route randomly */
	for (int r = 0; r < this->nb_vehicles; r++)
	{
		if (CL.empty())
		{
			route[r] = {0, 0};
			continue;
		}
		route[r] = {0, CL.back(), 0};
		load[r] = data->client[CL.back()].demand;
		CL.pop_back();
	}

	int cheapestInsertion = rand() % 2;

	bool forceFeasible = true;
	std::vector<int> routeAvailable(this->nb_vehicles, 1);
	int nbFullRoutes = 0;

	double gamma = (rand() % 35);
	gamma = gamma / 20.0;

	/** Pick a client and insert it into the best position possible until CL is empty */
	int r = 0;
	while (!CL.empty())
	{
		bool foundClient = false;
		double bestInsertionCost = std::numeric_limits<double>::infinity();
		auto bestNode = CL.begin();
		int bestPosition;

		if (routeAvailable[r] == 1)
		{
			for (int pos = 0; pos < route[r].size() - 1; pos++)
			{
				for (auto p = CL.begin(); p != CL.end(); ++p)
				{
					double insertionCost;

					if (cheapestInsertion)
					{
						insertionCost = data->time_cost[route[r][pos]][*p] + data->time_cost[*p][route[r][pos + 1]] - data->time_cost[route[r][pos]][route[r][pos + 1]] -
							gamma * (data->time_cost[0][*p] + data->time_cost[*p][0]);
					}
					else
					{
						insertionCost = data->time_cost[route[r][pos]][*p];
					}

					if (insertionCost < bestInsertionCost)
					{
						double totalDemand = load[r] + data->client[*p].demand;
						std::vector<double> sorted_deviations(route[r].size() - 1);
						for (int i = 1; i < route[r].size() - 1; i++)
						{
							sorted_deviations[i - 1] = data->client[route[r][i]].demand_deviation;
							// totalDemand += data->client[route[r][i]].demand;
						}
						sorted_deviations.back() = data->client[route[r][pos]].demand_deviation;

						if (data->budget_demand > 0)
						{
							if (route[r].size() - 2 >= data->budget_demand)
								std::partial_sort(sorted_deviations.begin(), sorted_deviations.begin(), sorted_deviations.end());
							for (int i = 0; i < sorted_deviations.size() && i < data->budget_demand; i++)
								totalDemand += sorted_deviations[i];
						}

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
				route[r].insert(route[r].begin() + bestPosition, *bestNode);
				load[r] += data->client[(*bestNode)].demand;
				CL.erase(bestNode);
			}
			else
			{
				routeAvailable[r] = 0;
				nbFullRoutes++;
				// forceFeasible = false;
			}
		}
		if (r == nb_vehicles - 1)
		{
			r = -1;
		}
		r++;
	}
}


void Solution::SetInterNeighbStatusAllTrue(int r1, int r2)
{
	for (int r = 0; r < route.size(); r++)
	{
		for (int n = 0; n < 6; n++)
		{
			inter_neighb_status[r1][r][n] = true;
			inter_neighb_status[r2][r][n] = true;
			inter_neighb_status[r][r1][n] = true;
			inter_neighb_status[r][r2][n] = true;
		}
	}
}

void Solution::SetInterNeighbStatusAllTrue(int r1)
{
	for (int r = 0; r < route.size(); r++)
	{
		for (int n = 0; n < 7; n++)
		{
			inter_neighb_status[r1][r][n] = true;
			inter_neighb_status[r][r1][n] = true;
		}
	}
}

void Solution::SetInterNeighbStatus(int r1, int r2, int neighb, bool status)
{
	inter_neighb_status[r1][r2][neighb] = status;
}

bool Solution::GetInterNeighbStatus(int r1, int r2, int neighb)
{
	return inter_neighb_status[r1][r2][neighb];
}

void Solution::SetIntraNeighbStatusAllTrue(int r)
{
	intra_neighb_status[r].assign(6, true);
}

void Solution::InitializeNeighbStatus()
{
	inter_neighb_status.assign(
			(int) route.size(), std::vector<std::vector<int>> (
				(int) route.size(), std::vector<int> (7, true)
				)	
			);

	intra_neighb_status.assign(
			(int) route.size(), std::vector<int> (6, true)			
			);
}

std::set<std::pair<int, int>> Solution::GetArcSet()
{
	std::set<std::pair<int, int>> arcset;
	for (int r = 0; r < route.size(); r++)
	{
		for (int i = 0; i < route[r].size() - 1; i++)
		{
			int client = route[r][i];
			int next = route[r][i+1];
			arcset.insert({client, next});
		}
	}

	return arcset;
}

int Solution::GetDistance(Solution *s)
{
	std::set<std::pair<int,int>> arcset1 = this->GetArcSet();
	std::set<std::pair<int,int>> arcset2 = s->GetArcSet();
	std::set<std::pair<int,int>> res;

	std::set_symmetric_difference(arcset1.begin(), arcset1.end(), arcset2.begin(), arcset2.end(), std::inserter(res, res.begin()));
	//	std::cout << res.size();
	return res.size();
}
