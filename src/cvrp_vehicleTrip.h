#ifndef CVRP_VEHICLE_TRIP
#define CVRP_VEHICLE_TRIP

#include <list>
#include "cvrp_idataModel.h"

namespace cvrp
{
class VehicleTrip
{
    public:
        VehicleTrip(const IDataModel& model);

        double cost() const { return m_cost; }
        int demandCovered() const { return m_demandCovered; }

        bool canAccommodate(int clientId) const;
        void addClientToTrip(int clientId);
        void optimiseCost();
        void reEvaluateDemandAndCost();
        bool isValidTrip() const;
        const std::list<int>& clientSeqConst() const { return m_clientSequence; }
        std::list<int>& clientSequence() { return m_clientSequence; }
	size_t getSeqSize() const { return m_clientSequence.size(); }
        std::string getTripStr() const;

    private:
        std::list<int> m_clientSequence;
        double m_cost;
        int m_demandCovered;
        const IDataModel& m_model;
};

}//cvrp namespace
#endif
