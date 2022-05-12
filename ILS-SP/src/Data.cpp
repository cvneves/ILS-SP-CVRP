#include "Data.h"

Data::Data(std::map<std::string, std::string> &params)
{
	if (params.find("-seed") != params.end())
	{
		seed = std::stoi(params["-seed"]);
		srand(seed);
	}
	else
	{
		seed = time(NULL);
		std::cout << "Seed: " << seed << std::endl;
		srand(seed);
	}

	GetName(params["-i"]);

	variant = params["-var"];
	if (variant == "rvrptw")
	{
		LoadInstanceVRPTW(params["-i"]);
		budget_demand = std::atoi(params["-bq"].c_str());
		budget_time = std::atoi(params["-bt"].c_str());
		std::string s = params["-aq"];
		std::replace(s.begin(), s.end(), ',', '.');
		// std::cout << s << "\n\n\n";
		alpha_demand = std::stod(s.c_str());
		s = params["-at"];
		std::replace(s.begin(), s.end(), ',', '.');
		alpha_time = std::stod(s.c_str());
		ComputeRVRPTWParameters();
	}
	if (variant == "rcvrp")
	{

		LoadInstanceCVRP(params["-i"]);

		std::string::size_type loc = name.find_last_of("-", name.size());

		std::string nVehiStr;

		if (loc != std::string::npos)
		{
			nVehiStr.append(name, loc + 2, name.size());
		}
		nb_vehicles = std::atoi(nVehiStr.c_str());

		for (int i = 0; i <= nb_clients; i++)
		{
			client[i].index = i;
		}

#ifdef ROBUST_DEMAND
		ComputeRCVRPParameters();
#endif
	}
	// CalculateProximity();
	// ComputeRandomCost();
	// PreProcessFsr();

	if (params.find("-ac") != params.end())
	{
		
		if (std::atoi(params["-ac"].c_str()) == 0)
		{
			use_acceptance_criterion = false;
		}
		else
			use_acceptance_criterion = true;
	}

	if (params.find("-ss") != params.end())
	{
		if (std::atoi(params["-ss"].c_str()) == 0)
			use_swap_star = false;
		else
			use_swap_star = true;
	}
	
	if (params.find("-ro") != params.end())
	{
		if (params["-ro"] == "1")
		{
			dist_type = 1;
		}
		else if (params["-ro"] == "0")
		{
			dist_type = 0;
		}
		else
		{
			dist_type = 2;
		}
	}
}

void Data::ComputeRandomCost()
{
	random_cost = std::vector<std::vector<double>>(nb_clients + 1, std::vector<double>(nb_clients + 1, 0));
	for (int i = 0; i < nb_clients + 1; i++)
	{
		for (int j = 0; j < nb_clients + 1; j++)
		{

			random_cost[i][j] = (double) rand() / (RAND_MAX + 1.0);
			// std::cout << random_cost[i][j] << " ";
			random_cost[i][j] *= 1;
			if (i == j)
				random_cost[i][j] = 0;
		}
		// std::cout << "\n";
	}
}

void Data::ComputeRCVRPParameters()
{
	budget_demand = floor((double)avg_route_len_factor * (double)(nb_clients) / nb_vehicles);

	min_cap = FindMinFeasibleCapacity(nb_vehicles);
	max_cap = FindMinFeasibleCapacity(nb_vehicles - 1);

	capacity = floor(cap_weight * max_cap + (1 - cap_weight) * min_cap + EPSILON);

	for (int i = 0; i < nb_clients + 1; i++)
	{
		client[i].demand_deviation = (client[i].demand * dem_dev_factor);
	}
}

int Data::CalcMinNbVehicles()
{
	int nbCustomers = nb_clients + 1;
	std::vector<int> sortedDemands(nbCustomers, 0);

	for (int i = 0; i < nbCustomers; i++)
	{
		sortedDemands[i] = client[i].demand;
	}

	std::sort(sortedDemands.begin(), sortedDemands.end());

	if (((1.0 + dem_dev_factor) * double(sortedDemands[nbCustomers - 1])) > capacity)
		return 999999999;

	int minVehicles = 1;
	double load = 0;
	int numDevs = 0;

	for (int i = nbCustomers - 1; i >= 0; i--)
	{
		double dem = double(sortedDemands[i]);

		if (numDevs < budget_demand)
			dem += dem_dev_factor * double(sortedDemands[i]);

		if ((load + dem) > capacity)
		{
			for (int j = i - 1; j >= 0; j--)
			{
				dem = double(sortedDemands[j]);
				if (numDevs < budget_demand)
					dem += dem_dev_factor * double(sortedDemands[j]);

				if ((load + dem) <= capacity)
				{
					int aux = sortedDemands[j];
					for (int k = j; k < i; k++)
						sortedDemands[k] = sortedDemands[k + 1];
					sortedDemands[i] = aux;
					break;
				}
			}

			if ((load + dem) > capacity)
			{
				minVehicles++;
				load = 0;
				numDevs = 0;
				dem = double(sortedDemands[i]);
				if (numDevs < budget_demand)
					dem += dem_dev_factor * double(sortedDemands[i]);
			}
		}

		load += dem;
		if (numDevs < budget_demand)
			numDevs++;
	}

	// std::cout << minVehicles << "\n";
	return minVehicles;
}

int Data::FindMinFeasibleCapacity(int nb_vehicles)
{
	int min_cap, max_cap;

	int sumDem = 0;
	int maxDem = 0;

	for (int i = 1; i < nb_clients + 1; i++)
	{
		sumDem += client[i].demand;
		if (client[i].demand > maxDem)
			maxDem = client[i].demand;
	}

	min_cap = (sumDem + nb_vehicles - 1) / nb_vehicles;
	max_cap = int(ceil(double(min_cap + maxDem) * (1.0 + dem_dev_factor)) + 1e-7);
	min_cap--;

	while ((max_cap - min_cap) > 1)
	{
		int testCap = (max_cap + min_cap) / 2;
		capacity = testCap;
		if (CalcMinNbVehicles() > nb_vehicles)
			min_cap = testCap;
		else
			max_cap = testCap;
	}

	return max_cap;
}

void Data::CalculateProximity()
{
	double gammaTW = 1.0;
	double gammaWT = 0.2;

	correlated_vertices = std::vector<std::vector<int>>(nb_clients + 1);

	std::vector<std::vector<std::pair<double, int>>> clientCorrelation(nb_clients + 1, std::vector<std::pair<double, int>>(nb_clients + 1));

	for (int i = 0; i <= nb_clients; i++)
	{
		for (int j = 0; j <= nb_clients; j++)
		{
			if (i == j)
				clientCorrelation[i][j].first = std::numeric_limits<double>::infinity();
			else
				clientCorrelation[i][j].first = time_cost[j][i] + gammaWT * std::max(client[i].ready_time - client[j].serv_time - time_cost[j][i] - client[j].due_time, 0.0) + gammaTW * std::max(client[j].ready_time + client[j].serv_time + time_cost[j][i] - client[i].due_time, 0.0);
			// clientCorrelation[i][j].first = time_cost[i][j] + gammaWT * std::max(client[j].ready_time - client[i].serv_time - time_cost[i][j] - client[i].due_time, 0.0) + gammaTW * std::max(client[i].ready_time + client[i].serv_time + time_cost[i][j] - client[j].due_time, 0.0);
			clientCorrelation[i][j].second = j;
		}
		std::sort(clientCorrelation[i].begin(), clientCorrelation[i].end());
		clientCorrelation[i].pop_back();
		for (int j = 0, k = 0; k < std::min(granular_size, nb_clients - 1); j++)
		{
			if (clientCorrelation[i][j].second == 0)
				continue;
			correlated_vertices[i].push_back(clientCorrelation[i][j].second);
			k++;
		}
		// cout << i << ": ";

		// for (int j = 0; j < clientCorrelation[i].size(); j++)
		// 	cout << clientCorrelation[i][j].second << " ";
		// cout << endl;
	}

	// for (int i = 0; i < correlated_vertices.size(); i++)
	// {
	// 	for (int j = 0; j < correlated_vertices[i].size(); j++)
	// 	{
	// 		std::cout << correlated_vertices[i][j] << " ";
	// 	}
	// 	std::cout << std::endl;
	// }
}

void Data::PreProcessFsr()
{
	ins_before.resize(nb_clients + 1, std::vector<int> (nb_clients + 1, 1));
	
	for (int i = 1; i <= nb_clients; i++)
		for (int j = 1; j <= nb_clients; j++)
		{
			if (client[i].ready_time + client[i].serv_time + time_cost[i][j] > client[j].due_time + EPSILON)
			{
				ins_before[i][j] = 0;
			}
		}
}

void Data::ComputeRVRPTWParameters()
{
	time_deviation = std::vector<std::vector<double>>(nb_clients + 1);

	for (int i = 0; i < nb_clients + 1; i++)
	{
		time_deviation[i] = (std::vector<double>(nb_clients + 1));
		for (int j = 0; j < nb_clients + 1; j++)
		{
			time_deviation[i][j] = 0.1 * floor(alpha_time * 10 * time_cost[i][j]);
		}
	}

	for (int i = 1; i < nb_clients + 1; i++)
	{
		client[i].demand_deviation += floor(alpha_demand * client[i].demand);
	}
}

void Data::ComputeDistanceMatrix()
{
	if (variant == "rvrptw" || dist_type == 2)
		return;
	
	time_cost = std::vector<std::vector<double>>(nb_clients + 1, std::vector<double>(nb_clients + 1));

	// Calcular Matriz Distancia (Euclidiana)
	for (int i = 0; i < nb_clients + 1; i++)
	{
		for (int j = 0; j < nb_clients + 1; j++)
		{
			//time_cost[i][j] = floor ( CalcDistEuc ( client[i].x, client[j].y, i, j ) + 0.5 );
			if (dist_type == 0)
				time_cost[i][j] = std::sqrt(std::pow(client[i].x - client[j].x, 2) + std::pow(client[i].y - client[j].y, 2));
			else if (dist_type == 1)
				time_cost[i][j] = std::floor(std::sqrt(std::pow(client[i].x - client[j].x, 2) + std::pow(client[i].y - client[j].y, 2)) + 0.5);
		}
	}
}

void Data::ComputeAuxVectors()
{
	demand = std::vector<double>();
	demand_deviation = std::vector<double>();

	for (int i = 0; i <= nb_clients; i++)
	{
		demand.push_back(client[i].demand);
		demand_deviation.push_back(client[i].demand_deviation);
	}
}

// Funcao LoadInstance VRPTW
void Data::LoadInstanceVRPTW(std::string &fileName)
{
	char instancia[100];
	char tempString[100];

	std::ifstream Instancia(fileName, std::ios::in);

	if (!Instancia)
	{
		std::cerr << "Arquivo nao pode ser aberto" << std::endl;
	}

	Instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia;

	Instancia >> nb_clients;
	
	//std::cout << nb_clients;
	//std::cout << "\n\n\n";

	while (atoi(instancia) == nb_clients)
	{
		Instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia;
		Instancia >> instancia;
		nb_clients++;
	}
	nb_clients--;

	Instancia.clear();
	Instancia.seekg(std::ios::beg);

	Instancia >> instancia >> instancia >> instancia >> instancia;
	Instancia >> nb_vehicles;
	Instancia >> capacity;

	Instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia >> instancia;

	client = std::vector<Vertex>(nb_clients + 1);

	for (int i = 0; i < nb_clients + 1; i++)
	{
		Instancia >> instancia >> client[i].x >> client[i].y >> client[i].demand >> client[i].ready_time >> client[i].due_time >> client[i].serv_time;
	}

	time_cost = time_deviation = std::vector<std::vector<double>>(nb_clients + 1, std::vector<double>(nb_clients + 1));

	for (int i = 0; i < nb_clients + 1; i++)
	{
		for (int j = 0; j < nb_clients + 1; j++)
		{
			//         matrizDistancia[i][j] = floor ( sqrt ( pow ( coord_x[i] - coord_x[j], 2 ) + pow ( coord_y[i] - coord_y[j], 2 ) ) + 0.5 );
			//      	 matrizDistancia[i][j] = sqrt ( pow ( coord_x[i] - coord_x[j], 2 ) + pow ( coord_y[i] - coord_y[j], 2 ) );

			//truncar na primeira casa decimal
			time_cost[i][j] = 10 * (sqrt(pow(client[i].x - client[j].x, 2) + pow(client[i].y - client[j].y, 2)));
			time_cost[i][j] = floor(time_cost[i][j]);
			time_cost[i][j] = time_cost[i][j] / 10;
		}
	}

	Instancia.close();

	// cout << "fim leitor VRPTW" << endl;
}

void Data::GetName(std::string &instanceName)
{

	std::string::size_type loc = instanceName.find_last_of(".", instanceName.size());
	std::string::size_type loc2 = instanceName.find_last_of("/", instanceName.size());

	if (loc != std::string::npos)
	{
		name.append(instanceName, loc2 + 1, loc - loc2 - 1);
	}
	else
	{
		name.append(instanceName, loc2 + 1, instanceName.size());
	}
}

void Data::LoadInstanceCVRP(std::string &fileName)
{
	// Open File

	std::ifstream Instance(fileName, std::ios::in);

	if (!Instance)
	{
		std::cerr << "File not found" << std::endl;
	}

	std::string temp, ewt;
	Instance >> temp;

	if (temp.compare("NAME") == 0)
	{
		Instance >> temp;

		while (temp.compare("DIMENSION:") != 0 && temp.compare("DIMENSION") != 0)
		{
			Instance >> temp;
		}

		if (temp.compare("DIMENSION") == 0)
			Instance >> temp;

		Instance >> nb_clients;
		nb_clients--;

		while (temp.compare("EDGE_WEIGHT_TYPE:") != 0 && temp.compare("EDGE_WEIGHT_TYPE") != 0)
		{
			Instance >> temp;
		}
		if (temp.compare("EDGE_WEIGHT_TYPE") == 0)
			Instance >> temp;

		Instance >> ewt;

		// Creating the Distance Matrix (depot = n)

		// Creating the demand arrays

		client = std::vector<Vertex>(nb_clients + 1);

		if (ewt == "EXPLICIT")
		{

			while (temp.compare("EDGE_WEIGHT_FORMAT:") != 0 && temp.compare("EDGE_WEIGHT_FORMAT") != 0)
			{
				Instance >> temp;
			}

			std::string ewf;
			if (temp.compare("EDGE_WEIGHT_FORMAT") == 0)
				Instance >> temp;
			Instance >> ewf;

			while (temp.compare("CAPACITY:") != 0 && temp.compare("CAPACITY") != 0)
				Instance >> temp;

			if (temp.compare("CAPACITY") == 0)
				Instance >> temp;

			Instance >> capacity;

			while (temp.compare("EDGE_WEIGHT_SECTION") != 0)
			{
				Instance >> temp;
			}

			time_cost = std::vector<std::vector<double>>(nb_clients + 1, std::vector<double>(nb_clients + 1));

			// Preencher Matriz Distancia
			for (int i = 1; i < nb_clients + 1; i++)
			{
				for (int j = 0; j < i; j++)
				{
					Instance >> time_cost[i][j];
					time_cost[j][i] = time_cost[i][j];
				}
			}

			for (int i = 0; i < nb_clients + 1; i++)
			{
				time_cost[i][i] = 0;
			}

			while (temp.compare("NODE_COORD_SECTION") != 0)
			{
				Instance >> temp;
			}
			// ler coordenadas
			for (int i = 0; i < nb_clients + 1; i++)
			{
				Instance >> temp;
				Instance >> client[i].x;
				Instance >> client[i].y;
			}
		}

		else if (ewt == "EUC_2D")
		{
			Instance >> temp;

			if (temp.compare("CAPACITY") == 0)
				Instance >> temp;

			Instance >> capacity;

			while (temp.compare("NODE_COORD_SECTION") != 0)
			{
				Instance >> temp;
			}
			// ler coordenadas
			for (int i = 0; i < nb_clients + 1; i++)
			{
				Instance >> temp;
				Instance >> client[i].x;
				Instance >> client[i].y;
			}

			//// Calcular Matriz Distancia (Euclidiana)
			//for (int i = 0; i < nb_clients + 1; i++)
			//{
			//	for (int j = 0; j < nb_clients + 1; j++)
			//	{
			//		//dist[i][j] = floor ( CalcDistEuc ( x, y, i, j ) + 0.5 );
			//		time_cost[i][j] = std::floor(std::sqrt(std::pow(client[i].x - client[j].x, 2) + std::pow(client[i].y - client[j].y, 2)) + 0.5);
			//	}
			//}
		}

		//Read demand section

		Instance >> temp;

		for (int i = 0; i < nb_clients + 1; i++)
		{
			Instance >> temp;
			Instance >> client[i].demand;
		}

		Instance.close();
	}

	else
	{
		Instance >> capacity;
		Instance >> nb_clients;
		nb_clients--;

		// Creating the Distance Matrix (depot = n)

		client = std::vector<Vertex>(nb_clients + 1);

		// time_cost = std::vector<std::vector<double>>(nb_clients + 1, std::vector<double>(nb_clients + 1));

		for (int i = 0; i < nb_clients + 1; i++)
		{
			for (int j = 0; j < nb_clients + 1; j++)
			{
				Instance >> time_cost[i][j];
			}
		}

		for (int i = 1; i < N + 1; i++)
		{
			Instance >> client[i].demand;
		}
	}

	// Close file

	Instance.close();
}

// // Funcao LoadInstance PRP
// void Data::LoadInstancePRP(int argc, char **argv)
// {

// 	char *instancia;
// 	instancia = argv[1];
// 	char tempString[100];

// 	// Abre arquivo

// 	ifstream Instancia(argv[1], ios::in);

// 	if (!Instancia)
// 	{
// 		cerr << "Arquivo nao pode ser aberto" << endl;
// 	}

// 	// Obtem dados do arquivo

// 	Instancia >> N;		   // Obtem o número de nós
// 	Instancia >> VCW;	   // Obtem o peso do veículo
// 	Instancia >> C;		   // Obtem a capacidade do veículo
// 	Instancia >> minSpeed; // Obtem a velocidade mínima
// 	Instancia >> maxSpeed; // Obtem o velocidade máxima

// 	// Cria Matriz de Adjacencia

// 	matrizDistancia = new double *[N + 1]; // Matriz de Distancia

// 	for (int i = 0; i < N + 1; i++)
// 	{
// 		matrizDistancia[i] = new double[N + 1];
// 	}

// 	// Preenche a Matriz de Adjacencia a partir dos dados do arquivo

// 	for (int i = 0; i < N + 1; i++)
// 	{
// 		for (int j = 0; j < N + 1; j++)
// 		{
// 			Instancia >> matrizDistancia[i][j];
// 		}
// 	}

// 	// Cria Vetor Delivery, Vetor Pickup, Vetor Ready Time, Vetor Due Time e Vetor Service Time

// 	delivery = new double[N + 1];  // Vetor Delivery
// 	pickup = new double[N + 1];	   // Vetor Pickup
// 	ready_time = new double[N + 1]; // Vetor Ready Time
// 	due_time = new double[N + 1];   // Vetor Due Time
// 	servTime = new double[N + 1];  // Vetor Service Time

// 	// Preenche os Vetores delivery e pickup a partir dos dados do arquivo

// 	for (int i = 0; i < N + 1; i++)
// 	{
// 		Instancia >> instancia >> instancia >> delivery[i] >> ready_time[i] >> due_time[i] >> servTime[i];
// 		pickup[i] = 0;
// 	}

// 	rSizeMax = N;

// 	// Fecha o arquivo

// 	Instancia.close();

// 	cout << "fim leitor PRP" << endl;
// }
