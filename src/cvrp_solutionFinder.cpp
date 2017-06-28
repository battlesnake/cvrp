#include "cvrp_solutionFinder.h"
#include "cvrp_util.h"
#include <algorithm>
#include <random>

namespace cvrp
{

SolutionFinder::SolutionFinder(const IDataModel& model) : m_model(model), m_dnaSequence(model.getClients())
{
}

SolutionModel SolutionFinder::getNaiveSolution() const
{
	SolutionModel solution;
	solution.chromosomes().push_back(VehicleTrip(m_model));
	auto genome = m_dnaSequence;
	std::shuffle(genome.begin(), genome.end(), std::mt19937(std::random_device{}()));
	for (const auto& gene : genome)
	{
		bool clientOnTrip = false;
		for (auto& chromosome : solution.chromosomes())
		{
			if(chromosome.canAccommodate(gene))
			{
				chromosome.addClientToTrip(gene);
				clientOnTrip = true;
				break;
			}
		}
		if (!clientOnTrip)
		{
			solution.chromosomes().push_back(VehicleTrip(m_model));
			solution.chromosomes().back().addClientToTrip(gene);
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
	constexpr unsigned max_generations = 1000;
	constexpr unsigned mutations_per_generation = 1000;
	constexpr unsigned max_contiguous_null_generations = 10;
	constexpr unsigned initial_population = 1000;
	constexpr unsigned max_population = 1000;

	const bool progress = !getenv("HIDE_PROGRESS");
	const bool benching = getenv("BENCH");

	std::vector<SolutionModel> population;
	double leastCost = -std::numeric_limits<double>::infinity();
	unsigned null_generations = 0;
	unsigned end_count = 0;

	const auto increase_population = [&] () {
		population.emplace_back(getNaiveSolution());
		auto cost = population.back().getCost();
		if (cost < leastCost || std::isinf(leastCost)) {
			leastCost = cost;
		}
	};

	const auto& getBest = [&] () {
		SolutionModel *sol = nullptr;
		for (auto& p : population)
		{
			if (sol == nullptr || p.getCost() < sol->getCost())
			{
				sol = &p;
			}
		}
		return *sol;
	};

	/* Initial population */
	for (unsigned i = 0; i < initial_population; ++i)
	{
		increase_population();
	}

	printf("max_generations=%u, mutations_per_generation=%u, max_contiguous_null_generations=%u\n", max_generations, mutations_per_generation, max_contiguous_null_generations);

	for (unsigned generation_num = 0; generation_num < max_generations; ++generation_num)
	{
		std::vector<SolutionModel> generation;
		bool foundBetterGeneration = false;
		if (progress)
		{
			fprintf(stderr, "\rpopulation=%zu, round=%u/%u (%.1f%%), score=%.1f, null rounds=%u            ", population.size(), generation_num, max_generations, (generation_num * 100.0 / max_generations), leastCost, null_generations);
		}
		const auto costThreshold = getBest().getCost();
		#pragma omp parallel for
		for (size_t n = 0; n < population.size(); ++n)
		{
			const auto& oldSol = population[n];
			for (unsigned mutation = 0; mutation < mutations_per_generation; mutation++)
			{
				const auto newSol = make_crossover(oldSol);
				const auto currSolCost = newSol.getCost();
				if (currSolCost < costThreshold && newSol.isValid(m_model.numberOfClients()))
				#pragma omp critical
				{
					if (currSolCost < leastCost)
					{
						leastCost = currSolCost;
						foundBetterGeneration = true;
					}
					generation.push_back(newSol);
				}
			}
		}
		if (foundBetterGeneration)
		{
			if (generation.size() <= max_population)
			{
				population = std::move(generation);
			} else {
				auto end_it = generation.begin() + max_population;
				std::partial_sort(generation.begin(), end_it, generation.end(),
					[] (const auto& a, const auto& b) { return a.getCost() < b.getCost(); });
				population.clear();
				std::move(generation.begin(), end_it, std::back_inserter(population));
			}
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

	return getBest();
}

}//cvrp namespace
