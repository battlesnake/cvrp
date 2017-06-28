#include "cvrp_solutionFinder.h"
#include "cvrp_util.h"
#include <algorithm>
#include <random>
#include <omp.h>

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

	static thread_local std::mt19937 gen{std::random_device{}()};
	std::uniform_int_distribution<int> uniform(1, chromosomes.size() - 1);

	int crossoverSubject1 = uniform(gen);
	int crossoverSubject2 = uniform(gen);
	while (crossoverSubject2 == crossoverSubject1)
	{
		crossoverSubject2 = uniform(gen);
	}

	int smallChromosomeSize = chromosomes[crossoverSubject1].getSeqSize() < chromosomes[crossoverSubject2].getSeqSize()
			? chromosomes[crossoverSubject1].getSeqSize()
			: chromosomes[crossoverSubject2].getSeqSize();

	int crossoverPoint = std::uniform_int_distribution<int>(1, smallChromosomeSize - 1)(gen);

	auto& subject1 = chromosomes[crossoverSubject1].clientSequence();
	auto& subject2 = chromosomes[crossoverSubject2].clientSequence();

	if (uniform(gen) & 1)
	{
		Util::splitAndCascade(subject1, subject2, crossoverPoint);
	}
	else
	{
		Util::splitAndFlipCascade(subject1, subject2, crossoverPoint);
	}

	for (auto& chromosome : chromosomes)
	{
		chromosome.reEvaluateDemandAndCost();
	}
}

SolutionModel SolutionFinder::solutionWithEvolution() const
{
	constexpr unsigned max_generations = 1'000;
	constexpr unsigned mutations_per_generation = 100'000'000;
	constexpr unsigned max_contiguous_null_generations = 10;
	constexpr unsigned initial_population = 100'000;
	constexpr unsigned max_population = 100'000;

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
		std::vector<double> costs;
		costs.reserve(population.size());
		std::transform(population.cbegin(), population.cend(), std::back_inserter(costs),
			[] (const auto& p) { return p.getCost(); });
		const auto min_idx = std::min_element(costs.cbegin(), costs.cend()) - costs.cbegin();
		return population[min_idx];
	};

	/* Initial population */
	for (unsigned i = 0; i < initial_population; ++i)
	{
		increase_population();
	}

	printf("max_generations=%u, mutations_per_generation=%u, max_contiguous_null_generations=%u\ninitial_population=%u, max_population=%u\n", max_generations, mutations_per_generation, max_contiguous_null_generations, initial_population, max_population);

	for (unsigned generation_num = 0; generation_num < max_generations; ++generation_num)
	{
		std::vector<SolutionModel> generation;
		bool foundBetterGeneration = false;
		if (progress)
		{
			fprintf(stderr, "\rpopulation=%zu, round=%u/%u (%.1f%%), score=%.1f, null rounds=%u            ", population.size(), generation_num, max_generations, (generation_num * 100.0 / max_generations), leastCost, null_generations);
		}
		const auto costThreshold = getBest().getCost();
		const auto mutations_per_subject = mutations_per_generation / population.size();
		const bool parallel_outer = population.size() > (unsigned) omp_get_num_threads() * 20;
#pragma omp parallel for if(parallel_outer)
		for (size_t n = 0; n < population.size(); ++n)
		{
			const auto& oldSol = population[n];
#pragma omp parallel for if(!parallel_outer)
			for (unsigned mutation = 0; mutation < mutations_per_subject; mutation++)
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
			/* Cull to keep population limit if necessary */
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
			/* Take fittest half only */
			auto end_it = population.begin() + (population.size() / 2 + 1);
			std::partial_sort(population.begin(), end_it, population.end(),
				[] (const auto& a, const auto& b) { return a.getCost() < b.getCost(); });
			population.erase(end_it, population.end());
		}
	}

	if (progress)
	{
		fprintf(stderr, "\n");
	}

	return getBest();
}

}//cvrp namespace
