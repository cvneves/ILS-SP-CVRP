#ifndef USERCALLBACK_H
#define USERCALLBACK_H

#include <ilcplex/ilocplex.h>
#include <limits>
#include "Data.h"
#include "Solution.h"
#include "Pool.h"

class MyCutCallback : public IloCplex::UserCutCallbackI 
{
	public:
		MyCutCallback(IloEnv env, const IloBoolVarArray &x, const IloBoolVarArray &y, Pool *my_pool, Data *data);
		IloCplex::CallbackI* duplicateCallback() const;
		void main();
		IloModel model;

	private:
		Data *data;
		Pool *my_pool;
    IloBoolVarArray x; // permanent paths
    IloBoolVarArray y; // temp paths
};

#endif
