#ifndef CVRP_UTIL
#define CVRP_UTIL

#include <list>

namespace cvrp
{
class Util
{
    public:
        static double distance(int x1, int y1, int x2, int y2);
        static int generateRandomNumberInRange(int min, int max);
        static void splitAndCascade(std::list<int>& first, std::list<int>& second, int splitpoint);
        static void splitAndFlipCascade(std::list<int>& first, std::list<int>& second, int splitPoint);
        static void spliceAndCascade(std::list<int>& first, std::list<int>& second, int start, int end);
};

}//cvrp namespace
#endif
