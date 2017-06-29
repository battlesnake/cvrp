#include "cvrp_vehicleTrip.h"
#include "cvrp_util.h"

#include <sstream>

namespace cvrp
{

VehicleTrip::VehicleTrip(const IDataModel& model) : m_model(&model)
{
    m_demandCovered = 0;
    m_cost = 0.0;
}

std::string VehicleTrip::getTripStr() const
{
    std::stringstream stream;
    stream << "x->";
    for (std::vector<int>::const_iterator ite = m_clientSequence.begin();
                ite != m_clientSequence.end(); ++ite)
    {
        stream << *ite << "->";
    }
    stream << "x ------- " << m_demandCovered;
    return stream.str();
}

bool VehicleTrip::canAccommodate(int clientId) const
{
    return ((m_model->getClientDemand(clientId) + m_demandCovered) <= m_model->vehicleCapacity());
}

bool VehicleTrip::isValidTrip() const
{
    return m_demandCovered <= m_model->vehicleCapacity();
}

void VehicleTrip::addClientToTrip(int clientId)
{
    m_clientSequence.push_back(clientId);
    m_demandCovered += m_model->getClientDemand(clientId);
}

void VehicleTrip::reEvaluateDemandAndCost()
{
    m_demandCovered = 0;
    for (auto i : m_clientSequence)
    {
        m_demandCovered += m_model->getClientDemand(i);
    }
    optimiseCost();
}

size_t VehicleTrip::calcHash() const
{
    size_t ret = 1;
    for (auto x : m_clientSequence)
    {
        ret++;
		if (!x)
		{
			x = 1;
		}
        ret *= x;
        ret ^= x*x;
    }
    return ret;
}

void VehicleTrip::optimiseCost()
{
    m_cost = 0;
    for (unsigned int i = 0; i < m_clientSequence.size(); i++)
    {
        unsigned int nearestClientIndex = i;
        double leastCost;
        if (i == 0)
        {
            leastCost = m_model->getClientDistanceFromDepot(m_clientSequence[i]);
        }
        else
        {
            leastCost = m_model->distanceBetweenClients(m_clientSequence[i-1], m_clientSequence[i]);
        }
        for (unsigned int j = i+1; j < m_clientSequence.size(); j++)
        {
            m_hash *= m_clientSequence[j];
            double currCost;
            if (i==0)
            {
                currCost = m_model->getClientDistanceFromDepot(m_clientSequence[j]);
            }
            else
            {
                currCost = m_model->distanceBetweenClients(m_clientSequence[i-1], m_clientSequence[j]);
            }

            if (currCost < leastCost)
            {
                leastCost = currCost;
                nearestClientIndex = j;
            }
        }

        if (nearestClientIndex != i)
        {
            std::swap(m_clientSequence[i], m_clientSequence[nearestClientIndex]);
        }
        m_cost += leastCost;
    }
    m_cost += m_model->getClientDistanceFromDepot(m_clientSequence.back());
    m_hash = hash();
}

}//cvrp namespace
