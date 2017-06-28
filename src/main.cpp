#include <iostream>
#include <sstream>
#include <fstream>
#include "cvrp_dataModel.h"
#include "cvrp_solutionModel.h"
#include "cvrp_solutionFinder.h"
#include "cvrp_util.h"

using namespace cvrp;

int main(int argc, char *argv[])
{
    std::ifstream dataFile(argv[1], std::ifstream::binary);
    std::stringstream jsonStream;
    jsonStream << dataFile.rdbuf();
    DataModel model(jsonStream);
    SolutionFinder solutionFinder(model);

    SolutionModel solution = solutionFinder.solutionWithEvolution();

    solution.printSolution();
    std::cout << "Total Cost: " << solution.getCost() << std::endl;
    std::cout << "##################################################################" << std::endl;

    return 0;
}
