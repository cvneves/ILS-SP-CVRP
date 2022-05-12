#ifndef MOVEMEMO_H
#define MOVEMEMO_H

#include "Solution.h"
#include <map>
#include <unordered_map>

// enum SMD_STATUS
// {
// 	UNVISITED = -1, UNIMPROVED = 0,	IMPROVED = 1
// };
// 
// struct SMD
// {
// 	int r1; 
// 	int r2;
// 	int r1_pos1;
// 	int r1_pos2;
// 	int r2_pos1;
// 	int r2_pos2;
// 
// 	SMD_STATUS status;
// 
// 	double cost;
// 	bool visited;
// };
// 
// class MoveMemo
// {
// 	std::map<std::pair<int,int>, SMD> inter_mov_memo;
// 	std::unordered_map<double, SMD> intra_mov_memo;
// 
// 	public:
// 		SMD* getSMD(Solution *s, int r1, int r2, int neighborhood);
// 		SMD* getSMD(Solution *s, int r, int neighborhood);
// 		
// 		void push(Solution *s, int r1, int r2, int neighborhood);
// 		void push(Solution *s, int r, int neighborhood);
// };

#endif
