#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <iomanip>
#include <map>

struct ArgParser
{
	std::map<std::string, std::string> params;

	ArgParser(int argc, char *argv[])
	{
//#ifndef DIMACS
//		if (argc < 2 || (argc - 1) % 2 != 0)
//		{
//			std::cout << std::endl
//					  << "*********************************** ILSSP-RVRPTW ********************************" << std::endl;
//			std::cout << "Basic usage (essential parameters): ./rvrp [-i instancePath] [-var variant]" << std::endl
//					  << std::endl;
//			std::cout << std::setw(22) << std::left << "Parameter:"
//					  << "Description:" << std::endl;
//			std::cout << std::setw(22) << std::left << "[-i instancePath]"
//					  << "Path to instance file." << std::endl;
//			std::cout << std::setw(22) << std::left << "[-var variant]"
//					  << "VRP variant (e.g. rvrp, rvrptw)." << std::endl;
//			std::cout << std::setw(22) << std::left << "[-seed seedNumber]"
//					  << "Seed. Set to time(NULL) if not specified." << std::endl;
//			std::cout << std::setw(22) << std::left << "[-sol solutionPath]"
//					  << "Output path to best solution." << std::endl;
//			std::cout << std::setw(22) << std::left << "[-dot plotSolPath]"
//					  << "Output path to plotted solution in .dot format." << std::endl;
//			std::cout << std::setw(22) << std::left << "[-bq Gamma^q]"
//					  << "Demand budget" << std::endl;
//			std::cout << std::setw(22) << std::left << "[-bt Gamma^t]"
//					  << "Time budget" << std::endl;
//			std::cout << std::setw(22) << std::left << "[-aq Alpha^q]"
//					  << "Demand deviation factor" << std::endl;
//			std::cout << std::setw(22) << std::left << "[-at Alpha^t]"
//					  << "Time deviation factor" << std::endl;
//			// std::cout << "[-i instancePath]" << std::setw(50) << std::right<< "path to the instance file" << std::endl;
//			std::cout << "*********************************************************************************" << std::endl
//					  << std::endl;
//			// std::cout << ""
//			exit(0);
//		}
//		else
//		{
//			for (int i = 1; i < argc; i = i + 2)
//				params[argv[i]] = argv[i + 1];
//		}
// #else
			params["-aq"] = "0.25";
			params["-bq"] = "5";
			params["-at"] = "0.1";
			params["-bt"] = "0";
			// params["-ac"] = "1";
			// params["-ss"] = "1";
			// params["-var"] = "rvrptw";
			params["-var"] = "rcvrp";
			params["-i"] = argv[1];
			// params["-seed"] = "0";
			// params["-ro"] = argv[2];
			// params["-ro"] = "1";
// #endif
	}
};

#endif
