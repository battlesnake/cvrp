#include <iostream>
#include <sstream>
#include <fstream>
#include "cvrp_dataModel.h"
#include "cvrp_solutionModel.h"
#include "cvrp_solutionFinder.h"

using namespace cvrp;

int main(int argc, char *argv[])
{
    std::ifstream dataFile(argv[1], std::ifstream::binary);
    std::stringstream jsonStream;
    jsonStream << dataFile.rdbuf();
    DataModel model(jsonStream);
    SolutionFinder solutionFinder(&model);

    SolutionModel solution;
    solutionFinder.getSolution(solution);

    if (solutionFinder.validateSolution(solution))
    {
        solution.printSolution();
        std::cout << "Total Cost: " << solutionFinder.solutionCost(solution) << std::endl;
    }
    
    return 0;
}