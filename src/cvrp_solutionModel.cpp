#include "cvrp_solutionModel.h"
#include <iostream>

namespace cvrp
{

void SolutionModel::printSolution()
{
	for (std::vector<VehicleTrip>::const_iterator it = m_solution.begin();
			it != m_solution.end(); ++it)
	{
		std::cout << it->getTripStr() << std::endl;
	}
}

double SolutionModel::getCost() const
{
	double totalCost = 0.0;
	for (const auto& chromosome : chromosomesConst())
	{
		totalCost += chromosome.cost();
	}
	return totalCost;
}

bool SolutionModel::isValid(int num_clients) const
{
	std::vector<int> check(num_clients + 1, false);
	for (const auto& chromosome : chromosomesConst())
	{
		if (!chromosome.isValidTrip())
		{
			return false;
		}
		for (const auto& clientSeq : chromosome.clientSeqConst())
		{
			if (check[clientSeq])
			{
				return false;
			}
			check[clientSeq] = true;
		}
	}
	for (unsigned int i = 1; i < check.size(); i++)
	{
		if (!check[i])
		{
			return false;
		}
	}
	return true;
}

}//cvrp namespace
