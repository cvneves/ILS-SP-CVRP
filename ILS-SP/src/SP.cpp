#include "SP.h"

void SPSolverCplex::Run()
{
	int nb_used_vehicles = 0;
	best_nb_vehi = best_solution.nb_vehicles;
	for (int i = 0; i < best_solution.route.size(); i++)
	{
		if (best_solution.route[i].size() > 2)
			nb_used_vehicles++;
	}

	my_pool->ForceTemp(&best_solution);

	IloEnv env;
	IloModel model(env);
	IloObjective objective = IloMinimize(env);
	IloBoolVarArray x(env, my_pool->GetNbPermPath()); // temporary pool
	IloBoolVarArray y(env, my_pool->GetNbNonPermPath()); // permanent pool

	IloRangeArray exactly_one_route(env, data->nb_clients, 1.0, 1.0);

	IloExpr obj_expr(env);
	IloExpr nb_vehicle_expr(env);

	for (int i = 0; i < my_pool->GetNbPermPath(); i++)
	{
		auto route_info = my_pool->GetPermRoute(i);
		model.add(x[i]);
		obj_expr += route_info->cost * x[i];
		nb_vehicle_expr += x[i];
	}

	for (int i = 0; i < my_pool->GetNbNonPermPath(); i++)
	{
		auto route_info = my_pool->GetTempRoute(i);
		model.add(y[i]);
		obj_expr += route_info->cost * y[i];
		nb_vehicle_expr += y[i];
	}

	model.add(IloMinimize(env, obj_expr));
	model.add(nb_vehicle_expr == nb_used_vehicles);

	std::vector<std::vector<int>> perm_cust(data->nb_clients+1,
																					std::vector<int>());
	std::vector<std::vector<int>> temp_cust(data->nb_clients+1,
																					std::vector<int>());

	for (int i = 0; i < my_pool->GetNbPermPath(); i++)
	{
		auto route_info = my_pool->GetPermRoute(i);
		for (int j = 1; j < route_info->sequence.size() - 1; j++)
			perm_cust[route_info->sequence[j]].push_back(i);
	}

	for (int i = 0; i < my_pool->GetNbNonPermPath(); i++)
	{
		auto route_info = my_pool->GetTempRoute(i);
		for (int j = 1; j < route_info->sequence.size() - 1; j++)
			temp_cust[route_info->sequence[j]].push_back(i);
	}

	for (int i = 1; i < data->nb_clients + 1;	i++)
	{
		IloExpr sum(env);
		for (auto &j : perm_cust[i])
			sum += x[j];
		for (auto &j : temp_cust[i])
			sum += y[j];
		model.add(sum == 1);
	}


//	model.add(x);
//	model.add(y);

//	model.add(objective);
//	model.add(exactly_one_route);
//	model.add(atMostKVehicles);

//	IloRange atMostKVehicles(env, nb_used_vehicles, nb_used_vehicles);

	// for (int i = 0; i < my_pool->GetNbPermPath(); i++)
	// {
	// 	auto route = my_pool->GetPermRoute(i);

	// 	IloNumArray colCoeff(env, data->nb_clients);
	// 	for (int i = 0; i < data->nb_clients; i++)
	// 		colCoeff[i] = 0;
	// 	for (int j = 1; j < route->sequence.size() - 1; j++)
	// 		colCoeff[route->sequence[j] - 1] = 1;

	// 	// x.add(IloBoolVar(objective(route.second) + exactly_one_route(colCoeff)));

	// 	x.add(IloBoolVar(objective(route->cost) + exactly_one_route(colCoeff) + atMostKVehicles(1)));
	// }

	// for (int i = 0; i < my_pool->GetNbNonPermPath(); i++)
	// {
	// 	auto route = my_pool->GetTempRoute(i);

	// 	IloNumArray colCoeff(env, data->nb_clients);
	// 	for (int i = 0; i < data->nb_clients; i++)
	// 		colCoeff[i] = 0;
	// 	for (int j = 1; j < route->sequence.size() - 1; j++)
	// 		colCoeff[route->sequence[j] - 1] = 1;

	// 	// y.add(IloBoolVar(objective(route.second) + exactly_one_route(colCoeff)));

	// 	y.add(IloBoolVar(objective(route->cost) + exactly_one_route(colCoeff) + atMostKVehicles(1)));
	// }

	IloCplex cplex(model);
	cplex.setParam(IloCplex::TiLim, data->max_sp_time);
	// cplex.setParam(IloCplex::CutUp, best_solution.cost);
	cplex.setParam(IloCplex::Threads, 1);
	//cplex.setParam(IloCplex::Param::Preprocessing::Linear, 0);
	cplex.setParam(IloCplex::Param::MIP::Interval, 100);
	// cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1);
	//cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);

	MyIncumbentCallback *inc_cbk = new (env) MyIncumbentCallback(env, best_solution.cost, x, y, &best_solution, my_pool, data);
	MyBranchCallback *branch_cbk = new (env) MyBranchCallback(env, inc_cbk);
	MyCutCallback *cut_cbk = new (env) MyCutCallback(env, x, y, my_pool, data);

	cplex.use(inc_cbk);
	cplex.use(branch_cbk);
	// cplex.use(cut_cbk);

	// Adding MIP Start
	IloNumVarArray start_var(env);
	IloNumArray start_val(env);

	for (int r = 0, i = y.getSize() - 1; r < best_solution.route.size(); r++)
	{
		if (best_solution.route[r].size() <= 2)
			continue;

		start_var.add(y[i]);
		start_val.add(1.0);
		i--;
	}

	cplex.addMIPStart(start_var, start_val);

	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());

	auto before = std::chrono::system_clock::now();

	try
	{
		cplex.solve();
	}
	catch (IloException &e)
	{
		std::cerr << "Exception raised by cplex: " << e;
		exit(EXIT_FAILURE);
	}

	auto after = std::chrono::system_clock::now();
	std::chrono::duration<double> solve_time = after - before;

	std::cout << "SP time: " << solve_time.count() << std::endl;

	if (cplex.getCplexStatus() == IloCplex::Optimal 
			|| cplex.getCplexStatus() == IloCplex::OptimalTol
			|| cplex.getCplexStatus() == IloCplex::AbortUser)
	{
		ti_lim_exceeded = false;
		if (cplex.getIncumbentNode() == 0)
			solved_at_root = true;
		else
			solved_at_root = false;
	}
	else
	{
		ti_lim_exceeded = true;
		solved_at_root = false;
	}

	if (solved_at_root)
	{
		data->tolerance += 0.0005;
	}

	std::cout << "Updating max gap from " << my_pool->gap_tol;

	if ((solve_time.count() <= 11))
	{
		my_pool->gap_tol += 0.1;
	}

	if (ti_lim_exceeded)
	{
		data->tolerance -= 0.0005;
		my_pool->gap_tol -= 0.1;
		my_pool->gap_tol = std::max(0.5, my_pool->gap_tol);
	}

	std::cout << " to " << my_pool->gap_tol << std::endl;

	// IloCplex cp(cut_cbk->model);
	// cp.exportModel((data->name + ".lp").c_str());

	delete inc_cbk;
	delete branch_cbk;
	delete cut_cbk;

	start_val.end();
	start_var.end();

	env.end();
}

MyIncumbentCallback::MyIncumbentCallback(const IloEnv env, double initial_cost, const IloBoolVarArray &x, const IloBoolVarArray &y, Solution *best_solution, Pool *my_pool, Data *data) : IloCplex::IncumbentCallbackI(env), initial_cost(initial_cost), x(x), y(y), best_solution(best_solution), my_pool(my_pool), data(data)
{
	incumbent_solution = Solution(data);
}

IloCplex::CallbackI *MyIncumbentCallback::duplicateCallback() const
{
	return new (getEnv()) MyIncumbentCallback(getEnv(), initial_cost, x, y, best_solution, my_pool, data);
}
void MyIncumbentCallback::main()
{
	aborted_sp = false;

	// loading the incumbent solution

	/** Converting MIP solution */
	IloNumArray x_values(getEnv(), my_pool->GetNbPermPath());
	IloNumArray y_values(getEnv(), my_pool->GetNbNonPermPath());
	getValues(x_values, x);
	getValues(y_values, y);

	incumbent_solution.route = std::vector<std::vector<int>>();
	incumbent_solution.route_cost = std::vector<double>();
	incumbent_solution.route_load = std::vector<double>();
	incumbent_solution.route_stats = std::vector<RouteStats>();
	incumbent_solution.nb_vehicles = 0;

	// int rIdx = 0;

	int nb_from_perm = 0;
	int nb_from_temp = 0;
	std::vector<double> sol_gap_list;
	std::vector<int> nb_it_ils_list;
	std::vector<int> nb_restart_list;
	std::vector<char> perm_temp_list;
	int inc_nb_vehi = 0;

	for (int r = 0; r < x_values.getSize(); r++)
	{
		if (x_values[r] > 0.9)
		{
			inc_nb_vehi++;
			nb_from_perm++;
			auto p = my_pool->GetPermRoute(r);
			incumbent_solution.route.push_back(p->sequence);
			incumbent_solution.route_cost.push_back(p->cost);
			incumbent_solution.route_load.push_back(p->load);
			incumbent_solution.route_stats.push_back(p->stats);
			incumbent_solution.nb_vehicles++;
			//sol_gap_list.push_back(p->solution_gap);
			//nb_it_ils_list.push_back(p->nb_it_ils);
			//nb_restart_list.push_back(p->nb_restart);
			//perm_temp_list.push_back('p');
			//// incumbent_solution.route[rIdx] = p.first;
			//// incumbent_solution.route_cost[rIdx] = p.second;
			//// rIdx++;
		}
	}

	for (int r = 0; r < y_values.getSize(); r++)
	{
		if (y_values[r] > 0.9)
		{
			inc_nb_vehi++;
			nb_from_temp++;
			auto p = my_pool->GetTempRoute(r);
			incumbent_solution.route.push_back(p->sequence);
			incumbent_solution.route_cost.push_back(p->cost);
			incumbent_solution.route_load.push_back(p->load);
			incumbent_solution.route_stats.push_back(p->stats);
			incumbent_solution.nb_vehicles++;
			//sol_gap_list.push_back(p->solution_gap);
			//nb_it_ils_list.push_back(p->nb_it_ils);
			//nb_restart_list.push_back(p->nb_restart);
			//perm_temp_list.push_back('t');
			//// incumbent_solution.route[rIdx] = p.first;
			//// incumbent_solution.route_cost[rIdx] = p.second;
			//// rIdx++;
		}
	}

	incumbent_solution.cost = getObjValue();

	// if (best_nb_vehi > 20)
	// {
	// 	cout << best_nb_vehi << endl;
	// 	exit(0);
	// }

	for (int r = incumbent_solution.route.size(); r < best_nb_vehi; r++)
	{
		incumbent_solution.route.push_back({0, 0});
		incumbent_solution.route_load.push_back(0);
		incumbent_solution.route_cost.push_back(0);
		incumbent_solution.route_stats.push_back(RouteStats());
	}
	
	incumbent_solution.InitializeNeighbStatus();
	//for (; inc_nb_vehi < 

#ifdef DIMACS
	std::cout << incumbent_solution << "\n";
#endif

//	if (incumbent_solution.cost < best_solution->cost)
	if (incumbent_solution.cost < initial_cost)
	{
		// std::cout << std::endl
		// 	<< "Temp routes: " << (double) nb_from_temp / (nb_from_temp + nb_from_perm) * 100 
		// 	<< "%" << std::endl;
		// std::cout << "Perm routes: " << (double) nb_from_perm / (nb_from_temp + nb_from_perm) * 100 
		// 	<< "%" << std::endl;

		// std::cout << "Route -- Associated Sol. Gap -- Nb iter ils" << std::endl;

		// for (int i = 0; i < incumbent_solution.route.size(); i++)
		// {
		// 	if (incumbent_solution.route[i].size() > 2)
		// 	{
		// 		printf("(R%d%c: %.2f/%0.2f, %d, %d) ", 
		// 				i + 1, 
		// 				(incumbent_solution.route_stats[i].is_permanent ? 'p' : 't'),
		// 				incumbent_solution.route_stats[i].solution_gap, 
		// 				100 * (incumbent_solution.route_stats[i].solution_cost 
		// 					- IlsRvnd::best_of_all_cost) / IlsRvnd::best_of_all_cost, 
		// 				incumbent_solution.route_stats[i].nb_improvs,
		// 				incumbent_solution.route_stats[i].restart_id);
		// 	}
		// }

		// std::cout << std::endl;
		// std::cout << std::endl;

		IlsRvnd ilsRvnd(data, my_pool, &incumbent_solution);
		ilsRvnd.Run(1, data->nb_clients + 5 * data->nb_vehicles, data->tolerance);

		if (ilsRvnd.best.cost < best_solution->cost)
		{
			*best_solution = (ilsRvnd.best);
			IlsRvnd::best_of_all_cost = best_solution->cost;
#ifdef DIMACS
			std::cout << *best_solution << "\n";
#endif
		}
	}
	// std::cout << *best_solution << "\n";

	double lower_bound = getBestObjValue();

	if (lower_bound > best_solution->cost)
	{
		puts("Aborting SP");
		abort();
		aborted_sp = true;
	}
}

MyBranchCallback::MyBranchCallback(IloEnv env, MyIncumbentCallback *inc_cbk) : IloCplex::BranchCallbackI(env),
	inc_cbk(inc_cbk)
{
}

IloCplex::CallbackI *MyBranchCallback::duplicateCallback() const
{
	return new (getEnv()) MyBranchCallback(getEnv(), inc_cbk);
}

void MyBranchCallback::main()
{
	aborted_sp = false;
	double lower_bound = getBestObjValue();

	if (lower_bound >= inc_cbk->best_solution->cost)
	{
		abort();
		aborted_sp = true;
	}

	if (std::ceil(getObjValue()) >= inc_cbk->best_solution->cost)
	{
		prune();
	}
}
