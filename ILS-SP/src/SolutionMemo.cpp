#include "SolutionMemo.h"

double SolutionMemo::GetHash(Solution *s)
{
	return 0;
}


int SolutionMemo::find(Solution *s, int iteration)
{
	// hashVal = GetHash(s);
	// auto it = memo.find(hashVal); 
	// if (it == memo.end() || it->second == iteration)
	// 	return -1;

	// return it->second;
}


bool SolutionMemo::push(Solution *s, int iteration)
{
	if(find(s, iteration) != -1)
		return false;

	// memo[hashVal] = iteration;
	return true;
}
