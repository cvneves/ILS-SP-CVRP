#ifndef CONSTRUCTION_H
#define CONSTRUCTION_H

#include "Data.h"
#include "Solution.h"

bool ParallelInsertion(Data *data, Solution *solution);
bool SequentialInsertion(Data *data, Solution *solution);
Solution GenerateInitialSolution(Data *data, int maxIter, int nb_vehicles);

#endif
