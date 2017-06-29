#include "cvrp_solutionFinder.h"
#include "cvrp_util.h"
#include <algorithm>
#include <random>
#include <set>
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
	constexpr unsigned mutations_per_generation = 1'000'000;
	constexpr unsigned max_contiguous_null_generations = 10;
	constexpr unsigned initial_population = 100'000;
	constexpr unsigned max_population = 100'000;

	const bool progress = !getenv("HIDE_PROGRESS");
	const bool benching = getenv("BENCH");

	struct CostedSolution
	{
		double cost;
		SolutionModel model;
		CostedSolution(SolutionModel&& model) :
			cost(model.getCost()),
			model(std::move(model))
			{ }
		CostedSolution(const SolutionModel& model) :
			cost(model.getCost()),
			model(model)
			{ }
		bool operator < (const CostedSolution& other) const
			{ return cost < other.cost; }
	};

	std::set<CostedSolution> population;
	double leastCost = -std::numeric_limits<double>::infinity();
	unsigned null_generations = 0;
	unsigned end_count = 0;

	const auto increase_population = [&] () {
		const auto it = population.emplace(getNaiveSolution());
		if (!it.second) {
			return;
		}
		auto cost = it.first->cost;
		if (cost < leastCost || std::isinf(leastCost)) {
			leastCost = cost;
		}
	};

	const auto& getBest = [&] () {
		const CostedSolution *best_ptr = nullptr;
		for (const auto& sm : population)
		{
			if (!best_ptr || sm.cost < best_ptr->cost)
			{
				best_ptr = &sm;
			}
		}
		return *best_ptr;
	};

	/* Initial population */
	for (unsigned i = 0; i < initial_population; ++i)
	{
		increase_population();
	}

	printf("max_generations=%u, mutations_per_generation=%u, max_contiguous_null_generations=%u\ninitial_population=%u, max_population=%u\n", max_generations, mutations_per_generation, max_contiguous_null_generations, initial_population, max_population);

	for (unsigned generation_num = 0; generation_num < max_generations; ++generation_num)
	{
		std::vector<CostedSolution> generation;
		bool foundBetterGeneration = false;
		if (progress)
		{
			fprintf(stderr, "\rpopulation=%zu, round=%u/%u (%.1f%%), score=%.1f, null rounds=%u            ", population.size(), generation_num, max_generations, (generation_num * 100.0 / max_generations), leastCost, null_generations);
		}
		const auto prev_best = getBest();
		const auto mutations_per_subject = mutations_per_generation / population.size();
		const bool parallel_outer = population.size() > (unsigned) omp_get_num_threads() * 20;
		/* Buffer in contiguous container for simple parallelisation */
		std::vector<const CostedSolution *> contiguous;
		for (const auto& subject : population)
		{
			contiguous.emplace_back(&subject);
		}
#pragma omp parallel for if(parallel_outer)
		for (size_t n = 0; n < contiguous.size(); ++n)
		{
			const auto& oldSol = *contiguous[n];
#pragma omp parallel for if(!parallel_outer)
			for (unsigned mutation = 0; mutation < mutations_per_subject; mutation++)
			{
				auto newSol = CostedSolution(make_crossover(oldSol.model));
				if (newSol < prev_best && newSol.model.isValid(m_model.numberOfClients()))
#pragma omp critical
				{
					if (newSol.cost < leastCost)
					{
						leastCost = newSol.cost;
						foundBetterGeneration = true;
					}
					generation.emplace_back(std::move(newSol));
				}
			}
		}
		if (foundBetterGeneration)
		{
			/* Cull to keep population limit if necessary */
			if (generation.size() > max_population)
			{
				auto end_it = generation.begin() + max_population;
				std::partial_sort(generation.begin(), end_it, generation.end());
			}
			population.clear();
			size_t count = 0;
			for (auto& x : generation)
			{
				if (count++ == max_population)
				{
					break;
				}
				population.emplace(std::move(x));
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

	return std::move(getBest().model);
}

}//cvrp namespace
