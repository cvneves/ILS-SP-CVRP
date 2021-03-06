# ILS-SP-CVRP
This repository contains a reimplementation of the ILS-SP matheuristic [[1](#1), [3](#3)] for solving the Capacitated Vehicle Routing Problem. The method relies on the Iterated Local Search (ILS) [[4](#4)] metaheuristic for constructing and improving solutions. During each iteration, routes belonging to the best solutions are feeded to a pool. A MIP solver is then employed to find the best solution by combining routes stored in the pool, which reduces to the Set Partitioning (SP) Problem.

In its current state, the code is aimed specifically at the canonical version of the VRP, which is constrained only by the vehicle capacities. Furthermore, in order to address the more recent benchmark set proposed in [[2](#2)], slight changes were made to the previously cited versions of the algorithm. These changes consist mainly in tweaks in the algorithm's parameters (e.g., number of iterations, route pool acceptance tolerance) and allowing an unrestricted number of vehicles in the solutions (as opposed to a fixed fleet). Instances from the benchmark set just mentioned, taken from [CVRPLIB](http://vrp.galgos.inf.puc-rio.br/index.php/en/), were also included in the `Instances` folder.

## Requirements

Both GCC 9.4.0 and [IBM CPLEX 20.1](https://www.ibm.com/products/ilog-cplex-optimization-studio) are required to compile the code. The latter may be obtained freely through an academic license. The makefile is configured to import the CPLEX libraries from their default installation folder (`/opt/ibm/ILOG/CPLEX_StudioXXX`). Nevertheless, importing them from a different folder, as well as using prior versions of the solver, should require minor changes in the makefile only.

## Usage 

Once the repository is cloned, compile the code by going to the `ILS-SP` directory and using the `make` command. Finally, enter the command
```sh
./ilssp.out ../Instances/instance_name.vrp
```
to run.


## References

<a id="1">[1]</a> Subramanian, A.; Uchoa, E.; Ochi, L. (2013).
A hybrid algorithm for a class of vehicle routing problems. Computers & Operations Research, 40(10), 2519-2531. ([DOI](https://doi.org/10.1016/j.cor.2013.01.013)).

<a id="2">[2]</a> Uchoa, E.; Pecin, D.; Pessoa, A.; Poggi de Arag??o, M.; Vidal, T.; Subramanian, A. (2017). New benchmark instances for the Capacitated Vehicle Routing Problem. European Journal of Operational Research, 257(3), 845-858. ([DOI](https://doi.org/10.1016/j.ejor.2016.08.012))

<a id="3">[3]</a> Subramanian, A.; Ochi, L.; Uchoa, E. (2012). Heuristic, Exact and Hybrid Approaches for Vehicle Routing Problems. Phd Thesis. ([Download](http://www.ic.uff.br/PosGraduacao/frontend-tesesdissertacoes/download.php?id=532.pdf&tipo=trabalho))

<a id="4">[4]</a> H. Louren??o; O. Martin; T. St??tzle. Iterated Local Search: Framework and Applications. Handbook of Metaheuristics. M. Gendreau; J-Y. Potvin (Eds.). Springer International Publishing (2019), pp. 129-168. ([DOI](https://doi.org/10.1007/978-3-319-91086-4_5))
