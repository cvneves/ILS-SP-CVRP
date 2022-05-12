#ifndef POOL_H
#define POOL_H

#include <unordered_set>
#include <unordered_map>
#include "Data.h"
#include "Solution.h"
#include "RouteStats.h"

class Pool
{
	public:
		Pool(Data *data) : data(data) {}
		virtual void Update(Solution *s, bool is_permanent) = 0; // adds paths in the solution to the pool
		virtual void ForceTemp(Solution *s) = 0;
		virtual void RemoveTempRoutes() = 0;
		virtual RouteInfo *GetPermRoute(int p) = 0;
		virtual RouteInfo *GetTempRoute(int p) = 0;
		virtual int GetNbPermPath() = 0;
		virtual int GetNbNonPermPath() = 0;
		virtual void GapClear(double best_sol_cost) = 0;
		double gap_tol = 0.5;
		double gap_step = 0.01;

	protected:
		Data *data;
		// std::unordered_map<double, tRoute> permPool;
		// std::unordered_map<double, tRoute> tempPool;
};

class SetPool : public Pool
{
public:
    SetPool(Data *data) : Pool(data) {}
    void Update(Solution *s, bool is_permanent);
		void ForceTemp(Solution *s);
    void RemoveTempRoutes();
    RouteInfo *GetPermRoute(int p);
    RouteInfo *GetTempRoute(int p);
    int GetNbPermPath();
    int GetNbNonPermPath();
		void GapClear(double best_sol_cost);

    std::vector<RouteInfo> permanent;
    std::vector<RouteInfo> temp;
    std::unordered_set<double> p_key;
    std::unordered_set<double> np_key;
    double GetHash(std::vector<int> &sequence);

private:
};

#endif
