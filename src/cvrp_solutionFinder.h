#ifndef CVRP_SOLUTION_FINDER
#define CVRP_SOLUTION_FINDER

#include "cvrp_idataModel.h"
#include "cvrp_solutionModel.h"

namespace cvrp
{
class SolutionFinder
{
    public:
        SolutionFinder(const IDataModel& model);

        SolutionModel getNaiveSolution() const;
        bool validateSolution(const SolutionModel& solution) const;
        SolutionModel solutionWithEvolution() const;

    private:
        const IDataModel& m_model;
        const std::vector<int> m_dnaSequence;

        void crossover(SolutionModel& solution) const;
        SolutionModel make_crossover(const SolutionModel& solution) const;
};

}//cvrp namespace
#endif
