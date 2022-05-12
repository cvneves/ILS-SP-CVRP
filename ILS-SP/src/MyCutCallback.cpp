#include "MyCutCallback.h"


MyCutCallback::MyCutCallback(IloEnv env, const IloBoolVarArray &x, const IloBoolVarArray &y, Pool *my_pool, Data *data) : IloCplex::UserCutCallbackI(env)
{

}

IloCplex::CallbackI* MyCutCallback::duplicateCallback() const
{
	return new (getEnv()) MyCutCallback(getEnv(), x, y, my_pool, data);
}

void MyCutCallback::main()
{
	std::cout << "Entering cutcallback...\n";
	model = getModel();	
	std::cout << getNrows() << "\n";

	//IloNumArray x_val(getEnv(), x.getSize());
	// IloNumArray y_val(getEnv(), y.getSize());

}
