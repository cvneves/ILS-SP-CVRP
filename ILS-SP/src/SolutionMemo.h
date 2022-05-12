#ifndef SOLUTION_MEMO_H
#define SOLUTION_MEMO_H
#include <unordered_map> 
#include <algorithm>
#include "Solution.h"

class SolutionMemo
{
	
	std::unordered_map<double, int> memo;
	double GetHash(Solution *s);

	public:
		int find(Solution *s, int iteration);
		bool push(Solution *s, int iteration); // true if the solution wasn't found before pushing
};

#endif
