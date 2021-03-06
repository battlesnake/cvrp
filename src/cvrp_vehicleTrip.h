#ifndef CVRP_VEHICLE_TRIP
#define CVRP_VEHICLE_TRIP

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
        const std::vector<int>& clientSeqConst() const { return m_clientSequence; }
        std::vector<int>& clientSequence() { return m_clientSequence; }
        size_t getSeqSize() const { return m_clientSequence.size(); }
        std::string getTripStr() const;

        bool operator == (const VehicleTrip& other) const
            { return m_clientSequence == other.m_clientSequence && m_model == other.m_model; }

        bool operator < (const VehicleTrip& other) const
            { return m_clientSequence < other.m_clientSequence && m_model == other.m_model; }

        size_t hash() const { return m_hash; }

    private:
        std::vector<int> m_clientSequence;
        double m_cost;
        int m_demandCovered;
        const IDataModel *m_model;
        size_t m_hash;
        size_t calcHash() const;
};

}//cvrp namespace
#endif
