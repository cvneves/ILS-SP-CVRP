# ILSSP-CVRP
This repository contains a reimplementation of the ILSSP matheuristic [[1](#1), [3](#3)], devised for solving the Capacitated Vehicle Routing Problem. The method relies on the Iterated Local Search metaheuristic for constructing and improving solutions. During each iteration, routes belonging to the best solutions are feeded to a pool. The problem of finding the best solution by using routes stored in the pool, which reduces to the Set Partitioning Problem, is solved by employing a MIP.

In its current state, the code is aimed specifically at the canonical version of the VRP, which considers only the vehicle capacities. Furthermore, slight modifications were made to the previously cited versions in order to address the more recent benchmark set proposed by [[2](#2)]. Some of these changes include tweaks in the algorithm's parameters (e.g., number of iterations, route pool acceptance tolerance) and allowing an unrestricted number of vehicles in the solutions (as opposed to a fixed fleet).

## Requirements and usage

Both GCC 9.4.0 and [IBM CPLEX 20.1](https://www.ibm.com/products/ilog-cplex-optimization-studio) are required to compile the code. The latter may be obtained freely through an academic license. The makefile is configured to import the CPLEX libraries from their default installation folder (`/opt/ibm/ILOG/CPLEX_StudioXXX`). Nevertheless, importing them from a different folder, as well as using prior versions of the solver, should require minor changes in the makefile only.

Once the repository is cloned, compile the code by going to the `ILSSP` directory and using the `make` command. Finally, enter the command
```sh
./rvrp.out ../Instances/instance_name.vrp
```
to run.


## References

<a id="1">[1]</a> Subramanian, A.; Uchoa, E.; Ochi, L. (2013).
A hybrid algorithm for a class of vehicle routing problems. Computers & Operations Research, 40(10), 2519-2531. ([DOI](https://doi.org/10.1016/j.cor.2013.01.013)).

<a id="2">[2]</a> Uchoa, E.; Pecin, D.; Pessoa, A.; Poggi de Aragão, M.; Vidal, T.; Subramanian, A. (2017). New benchmark instances for the Capacitated Vehicle Routing Problem. European Journal of Operational Research, 257(3), 845-858. ([DOI](https://doi.org/10.1016/j.ejor.2016.08.012))

<a id="3">[3]</a> Subramanian, A.; Ochi, L.; Uchoa, E. (2012). Heuristic, Exact and Hybrid Approaches for Vehicle Routing Problems. Phd Thesis. ([Download](http://www.ic.uff.br/PosGraduacao/frontend-tesesdissertacoes/download.php?id=532.pdf&tipo=trabalho))
