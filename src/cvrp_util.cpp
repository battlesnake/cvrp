#include "cvrp_util.h"

#include <algorithm>
#include <math.h>
#include <random>
#include <ctime>
#include <cstdlib>


namespace cvrp
{

double Util::distance(int x1, int y1, int x2, int y2)
{
    return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
}

int Util::generateRandomNumberInRange(int min, int max)
{
    static thread_local std::default_random_engine gen(std::random_device{}());
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(gen);
}

std::list<int>::iterator list_ptr(std::list<int>::iterator it, size_t idx)
{
    while (idx--) {
        ++it;
    }
    return it;
}

void Util::splitAndCascade(std::list<int>& first, std::list<int>& second, int splitPoint)
{
    /* Swaps tails of the lists */
    std::list<int> tmp;
    tmp.splice(tmp.begin(), first, list_ptr(first.begin(), splitPoint), first.end());
    first.splice(first.end(), second, list_ptr(second.begin(), splitPoint), second.end());
    second.splice(second.end(), tmp, tmp.begin(), tmp.end());
}

void Util::splitAndFlipCascade(std::list<int>& first, std::list<int>& second, int splitPoint)
{
    /* Swaps tail of first with head of second list */
    std::list<int> tmp;
    tmp.splice(tmp.begin(), first, list_ptr(first.begin(), splitPoint), first.end());
    first.splice(first.end(), second, second.begin(), list_ptr(second.begin(), splitPoint));
    second.splice(second.end(), tmp, tmp.begin(), tmp.end());
}

void Util::spliceAndCascade(std::list<int>& first, std::list<int>& second, int start, int end)
{
    /* Swaps subrange of each list */
    auto first_begin = list_ptr(first.begin(), start);
    auto second_begin = list_ptr(second.begin(), start);
    auto first_end = list_ptr(first_begin, end - start);
    auto second_end = list_ptr(second_begin, end - start);

    std::list<int> tmp;
    tmp.splice(tmp.begin(), first, first_begin, first_end);
    first.splice(first.end(), second, second_begin, second_end);
    second.splice(second.end(), tmp, tmp.begin(), tmp.end());
}

}//cvrp namespace
