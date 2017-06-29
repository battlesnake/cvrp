#ifndef CVRP_UTIL
#define CVRP_UTIL

#include <vector>
#include <random>

namespace cvrp
{
class Util
{
    public:
        static void seed_prngs();
        static std::mt19937& get_prng();
        static double distance(int x1, int y1, int x2, int y2);
        static void splitAndCascade(std::vector<int>& first, std::vector<int>& second, int splitpoint);
        static void splitAndFlipCascade(std::vector<int>& first, std::vector<int>& second, int splitPoint);
};

}//cvrp namespace
#endif
