#ifndef SP_H
#define SP_H

#include "RouteStats.h"
#include "Data.h"
#include "Solution.h"
#include "Pool.h"
#include "IlsRvnd.h"
#include <ilcplex/ilocplex.h>
#include "MyCutCallback.h"
#include <cstring>

static bool aborted_sp = false;

class SPSolver
{

public:
    SPSolver(Data *data, Pool *my_pool, Solution *solution) : data(data), my_pool(my_pool), best_solution(*solution) {}
    virtual void Run() = 0;
    Solution &getBestSolution() { return best_solution; }
    bool getTiLimExceeded() { return ti_lim_exceeded; }
    bool getSolvedAtRoot() { return solved_at_root; }
		void SetSolution(Solution *solution) { best_solution = *solution; }

protected:
    Data *data;
    Pool *my_pool;
    Solution best_solution;
    bool ti_lim_exceeded;
    bool solved_at_root;
};

class SPSolverCplex : public SPSolver
{

public:
    SPSolverCplex(Data *data, Pool *my_pool, Solution *solution) : SPSolver(data, my_pool, solution) {}
    void Run();
};

class MyIncumbentCallback : public IloCplex::IncumbentCallbackI
{
public:
    double initial_cost;
    MyIncumbentCallback(const IloEnv env, double initial_cost, const IloBoolVarArray &x, const IloBoolVarArray &y, Solution *best_solution, Pool *my_pool, Data *data);
    IloCplex::CallbackI *duplicateCallback() const;
    void main();

    IloBoolVarArray x; // permanent paths
    IloBoolVarArray y; // temp paths
    Solution *best_solution;
    Solution incumbent_solution;
    Pool *my_pool;
    Data *data;
};

class MyBranchCallback : public IloCplex::BranchCallbackI 
{
   public:
      MyBranchCallback(IloEnv env, MyIncumbentCallback* inc_cbk);
      IloCplex::CallbackI *duplicateCallback() const; 
      void main();
      MyIncumbentCallback* inc_cbk;
};

#endif
