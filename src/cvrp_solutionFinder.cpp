#include "cvrp_solutionFinder.h"
#include "cvrp_util.h"
#include <algorithm>
#include <set>
#include <omp.h>
#include <csignal>
#include <atomic>

namespace cvrp
{

SolutionFinder::SolutionFinder(const IDataModel& model) : m_model(model), m_dnaSequence(model.getClients())
{
}

SolutionModel SolutionFinder::getNaiveSolution(const std::vector<int>& genome) const
{
	SolutionModel solution;
	solution.chromosomes().push_back(VehicleTrip(m_model));
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

	std::uniform_int_distribution<int> uniform(1, chromosomes.size() - 1);

	auto gen = Util::get_prng();

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

std::atomic_bool sigend{false};

void sigend_handler(int)
{
	sigend = true;
}

SolutionModel SolutionFinder::solutionWithEvolution() const
{
	constexpr unsigned long max_generations = 100;
	constexpr unsigned long max_mutations_per_generation = 10'000'000'000;
	constexpr unsigned long max_contiguous_null_generations = 3;
	constexpr unsigned long initial_population = 100'000;
	constexpr unsigned long max_population = 10'000'000;
	constexpr unsigned long max_mutations_per_subject = 100'000;

	const bool progress = !getenv("HIDE_PROGRESS");
	const bool benching = getenv("BENCH");

	std::signal(SIGINT, sigend_handler);
	std::signal(SIGTERM, sigend_handler);

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
		bool operator == (const CostedSolution& other) const
			{ return cost == other.cost && model == other.model; }
		bool operator < ( const CostedSolution& other) const
			{ return cost < other.cost || model < other.model; }
	};

	struct CostedSolutionHash
	{
		size_t operator () (const CostedSolution& x) const
		{
			union {
				size_t s;
				double d;
			} t;
			t.d = x.cost;
			return t.s ^ x.model.hash();
		}
	};

	struct CostedSolutionEqual
	{
		bool operator () (const CostedSolution& x, const CostedSolution& y) const
			{ return x == y; }
	};

	struct CostedSolutionCompare
	{
		bool operator () (const CostedSolution& x, const CostedSolution& y) const
			{ return x < y; }
	};

	using ResultSet = std::set<CostedSolution, CostedSolutionCompare>;

	ResultSet population;
	unsigned null_generations = 0;

	Util::seed_prngs();

	/* Initial population */
	printf("Initialising %'lu random solutions\n", initial_population);
	#pragma omp parallel
	{
		ResultSet buf;
		auto& prng = Util::get_prng();
		auto genome = m_dnaSequence;
#pragma omp for
		for (unsigned long i = 0; i < initial_population; ++i)
		{
			std::shuffle(genome.begin(), genome.end(), prng);
			const auto it = buf.emplace(getNaiveSolution(genome));
			if (!it.second)
			{
				continue;
			}
		}
#pragma omp critical
		{
			population.merge(std::move(buf));
		}
	}

	printf("max_generations=%'lu, max_mutations_per_generation=%'lu, max_contiguous_null_generations=%'lu\ninitial_population=%'lu, max_population=%'lu\n", max_generations, max_mutations_per_generation, max_contiguous_null_generations, initial_population, max_population);

	for (unsigned long generation_num = 0; generation_num < max_generations; ++generation_num)
	{
		ResultSet generation;
		if (progress)
		{
			fprintf(stderr, "\rpopulation=%'zu, round=%'lu/%'lu (%.1f%%), score=%.1f, null rounds=%'u            ", population.size(), generation_num, max_generations, (generation_num * 100.0 / max_generations), population.begin()->cost, null_generations);
		}
		const auto threshold = (--population.end())->cost;
		const auto mutations_per_subject = std::min<size_t>(max_mutations_per_generation / population.size(), max_mutations_per_subject);
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
			for (unsigned long mutation = 0; mutation < mutations_per_subject; mutation++)
			{
				if (sigend)
				{
					continue;
				}
				auto newSol = CostedSolution(make_crossover(oldSol.model));
				if (newSol.cost < threshold && newSol.model.isValid(m_model.numberOfClients()))
#pragma omp critical
				{
					if (generation.empty() || newSol.cost < (--generation.end())->cost)
					{
						if (generation.size() == max_population)
						{
							generation.erase(--generation.end());
						}
						generation.emplace(std::move(newSol));
					}
				}
			}
		}
		if (!generation.empty() && generation.begin()->cost < population.begin()->cost)
		{
			null_generations = 0;
			population = std::move(generation);
			if (progress)
			{
				fprintf(stderr, "\n");
			}
		} else {
			null_generations++;
			if (null_generations == max_contiguous_null_generations && !benching) {
				break;
			}
		}
		if (sigend)
		{
			break;
		}
	}

	if (progress)
	{
		fprintf(stderr, "\n");
	}

	return std::move(population.begin()->model);
}

}//cvrp namespace
