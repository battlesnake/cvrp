#ifndef CVRP_DATA_MODEL
#define CVRP_DATA_MODEL

#include "cvrp_idataModel.h"
#include "json/json.h"

namespace cvrp
{
class DataModel : public IDataModel
{
    public:
        DataModel(std::stringstream& jsonData);

        double distanceBetweenClients(int client1Id, int client2Id) const;
        int vehicleCapacity() const { return m_vehicleCpacity; }
        const Coord& depot() const { return m_depot; }
        int getClientDemand(int clientId) const;
        const Coord& getClientLocation(int clientId) const;
        int numberOfClients() const;
	std::vector<int> getClients() const;
        double getClientDistanceFromDepot(int clientId) const;

    private:
        Clients m_clientDemands;
        int m_vehicleCpacity;
        Coord m_depot;
        void populateData(Json::Value& jsonObj);
	const Client& getClient(int clientId) const { return m_clientDemands.at(clientId); }
};

}//cvrp namespace
#endif
