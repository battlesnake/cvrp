#include "cvrp_solutionFinder.h"
#include "cvrp_util.h"

namespace cvrp
{

SolutionFinder::SolutionFinder(const IDataModel& model) : m_model(model), m_dnaSequence(model.getClients())
{
}

SolutionModel SolutionFinder::getNaiveSolution() const
{
	SolutionModel solution;
	solution.chromosomes().push_back(VehicleTrip(m_model));
	for (unsigned int i = 0; i < m_dnaSequence.size(); i++)
	{
		bool clientOnTrip = false;
		for (auto& chromosome : solution.chromosomes())
		{
			if(chromosome.canAccommodate(m_dnaSequence[i]))
			{
				chromosome.addClientToTrip(m_dnaSequence[i]);
				clientOnTrip = true;
				break;
			}
		}
		if (!clientOnTrip)
		{
			solution.chromosomes().push_back(VehicleTrip(m_model));
			solution.chromosomes().back().addClientToTrip(m_dnaSequence[i]);
		}
	}
	for (auto& chromosome : solution.chromosomes())
	{
		chromosome.optimiseCost();
	}
	return solution;
}

SolutionModel SolutionFinder::make_crossover(const SolutionModel& solution) const
{
	SolutionModel sm = solution;
	crossover(sm);
	return sm;
}

void SolutionFinder::crossover(SolutionModel& solution) const
{
	auto& chromosomes = solution.chromosomes();

	int crossoverSubject1 = Util::generateRandomNumberInRange(1,chromosomes.size()-1);
	int crossoverSubject2 = Util::generateRandomNumberInRange(1,chromosomes.size()-1);
	while (crossoverSubject2 == crossoverSubject1)
	{
		crossoverSubject2 = Util::generateRandomNumberInRange(1,chromosomes.size()-1);
	}

	int smallChromosomeSize = chromosomes[crossoverSubject1].getSeqSize() < chromosomes[crossoverSubject2].getSeqSize()
			? chromosomes[crossoverSubject1].getSeqSize()
			: chromosomes[crossoverSubject2].getSeqSize();

	int crossoverPoint = Util::generateRandomNumberInRange(1,smallChromosomeSize-1);

	auto& subject1 = chromosomes[crossoverSubject1].clientSequence();
	auto& subject2 = chromosomes[crossoverSubject2].clientSequence();

	int randomNum = Util::generateRandomNumberInRange(1,10);

	/*
	 * Might be faster to use linked lists splice/join rather than vector for
	 * this split/cascade part, TODO see what effect this would have on
	 * performance of the rest of the program.
	 */

	if (randomNum < 8)
	{
		Util::splitAndCascade(subject1, subject2, crossoverPoint);
	}
	else //if (randomNum < 9)
	{
		Util::splitAndFlipCascade(subject1, subject2, crossoverPoint);
	}
	// else
	// {
	//     int crossoverPoint2 = Util::generateRandomNumberInRange(1,smallChromosomeSize-2);
	//     while (crossoverPoint == crossoverPoint2)
	//     {
	//         crossoverPoint2 = Util::generateRandomNumberInRange(1,smallChromosomeSize-2);
	//     }
	//     if (crossoverPoint2 < crossoverPoint)
	//     {
	//         Util::spliceAndCascade(subject1, subject2, crossoverPoint2, crossoverPoint);
	//     }
	//     else
	//     {
	//         Util::spliceAndCascade(subject1, subject2, crossoverPoint, crossoverPoint2);
	//     }
	// }

	for (auto& chromosome : chromosomes)
	{
		chromosome.reEvaluateDemandAndCost();
	}
}

SolutionModel SolutionFinder::solutionWithEvolution() const
{
	constexpr unsigned max_generations = 100;
	constexpr unsigned mutations_per_generation = 10000;
	constexpr unsigned max_contiguous_null_generations = 10;

	const bool progress = !getenv("HIDE_PROGRESS");
	const bool benching = getenv("BENCH");

	std::vector<SolutionModel> solutions;
	solutions.emplace_back(getNaiveSolution());
	auto leastCost = solutions.back().getCost();
	unsigned null_generations = 0;
	unsigned end_count = 0;

	printf("max_generations=%u, mutations_per_generation=%u, max_contiguous_null_generations=%u\n", max_generations, mutations_per_generation, max_contiguous_null_generations);

	for (unsigned generation_num = 0; generation_num < max_generations; ++generation_num)
	{
		std::vector<SolutionModel> generation;
		bool foundBetterGeneration = false;
		if (progress)
		{
			fprintf(stderr, "\rsolutions=%zu, round=%u/%u (%.1f%%), score=%.1f, null rounds=%u            ", solutions.size(), generation_num, max_generations, (generation_num * 100.0 / max_generations), leastCost, null_generations);
		}
		for (const auto& oldSol : solutions)
		{
			const auto startingCost = oldSol.getCost();
			#pragma omp parallel for
			for (unsigned mutation = 0; mutation < mutations_per_generation; mutation++)
			{
				const auto newSol = make_crossover(oldSol);
				const auto currSolCost = newSol.getCost();
				if (currSolCost < startingCost && newSol.isValid(m_model.numberOfClients()))
				#pragma omp critical
				{
					if (currSolCost < leastCost)
					{
						leastCost = currSolCost;
						foundBetterGeneration = true;
						generation.push_back(newSol);
					}
				}
			}
		}
		if (foundBetterGeneration)
		{
			solutions = std::move(generation);
			if (progress)
			{
				fprintf(stderr, "\n");
			}
		} else {
			null_generations++;
			if (end_count++ == max_contiguous_null_generations && !benching) {
				break;
			}
		}
	}

	if (progress)
	{
		fprintf(stderr, "\n");
	}

	return solutions.back();
}

}//cvrp namespace
