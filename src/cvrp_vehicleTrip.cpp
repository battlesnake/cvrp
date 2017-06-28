#include "cvrp_vehicleTrip.h"
#include "cvrp_util.h"

#include <sstream>
#include <functional>

namespace cvrp
{

VehicleTrip::VehicleTrip(const IDataModel& model) : m_model(model)
{
    m_demandCovered = 0;
    m_cost = 0.0;
}

std::string VehicleTrip::getTripStr() const
{
    std::stringstream stream;
    stream << "x->";
    for (const auto& waypoint : m_clientSequence)
    {
        stream << waypoint << "->";
    }
    stream << "x ------- " << m_demandCovered;
    return stream.str();
}

bool VehicleTrip::canAccommodate(int clientId) const
{
    return ((m_model.getClientDemand(clientId) + m_demandCovered) <= m_model.vehicleCapacity());
}

bool VehicleTrip::isValidTrip() const
{
    return m_demandCovered <= m_model.vehicleCapacity();
}

void VehicleTrip::addClientToTrip(int clientId)
{
    m_clientSequence.push_back(clientId);
    m_demandCovered += m_model.getClientDemand(clientId);
}

void VehicleTrip::reEvaluateDemandAndCost()
{
    m_demandCovered = 0;
    for (auto i : m_clientSequence)
    {
        m_demandCovered += m_model.getClientDemand(i);
    }
    optimiseCost();
}

void VehicleTrip::optimiseCost()
{
    std::vector<int *> ordered;
    for (auto& x : m_clientSequence)
    {
        ordered.push_back(&x);
    }
    m_cost = 0;
    for (unsigned int i = 0; i < ordered.size(); i++)
    {
        unsigned int nearestClientIndex = i;
        double leastCost;
        if (i == 0)
        {
            leastCost = m_model.getClientDistanceFromDepot(*ordered[i]);
        }
        else
        {
            leastCost = m_model.distanceBetweenClients(*ordered[i-1], *ordered[i]);
        }
        for (unsigned int j = i+1; j < m_clientSequence.size(); j++)
        {
            double currCost;
            if (i==0)
            {
                currCost = m_model.getClientDistanceFromDepot(*ordered[j]);
            }
            else
            {
                currCost = m_model.distanceBetweenClients(*ordered[i-1], *ordered[j]);
            }

            if (currCost < leastCost)
            {
                leastCost = currCost;
                nearestClientIndex = j;
            }
        }

        if (nearestClientIndex != i)
        {
            std::swap(*ordered[i], *ordered[nearestClientIndex]);
        }
        m_cost += leastCost;
    }
    m_cost += m_model.getClientDistanceFromDepot(m_clientSequence.back());
}

}//cvrp namespace
