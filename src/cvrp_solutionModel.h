#ifndef CVRP_SOLUTION_MODEL
#define CVRP_SOLUTION_MODEL

#include <vector>
#include "cvrp_vehicleTrip.h"

namespace cvrp
{
class SolutionModel
{
    public:
        std::vector<VehicleTrip>& chromosomes() { return m_solution; }
        const std::vector<VehicleTrip>& chromosomesConst() const { return m_solution; }
        void printSolution();
        double getCost() const;
        bool isValid(int num_clients) const;

        bool operator == (const SolutionModel& other) const
            { return m_solution == other.m_solution; }

        bool operator < (const SolutionModel& other) const
            { return m_solution < other.m_solution; }

        size_t hash() const;

    private:
        std::vector<VehicleTrip> m_solution;
};

}//cvrp namespace
#endif
