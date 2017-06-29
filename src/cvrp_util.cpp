#include "cvrp_util.h"

#include <algorithm>
#include <cmath>

namespace cvrp
{

static thread_local std::mt19937 _prng(std::random_device{}());

std::mt19937 Util::get_prng()
{
    return _prng;
}

double Util::distance(int x1, int y1, int x2, int y2)
{
    long dx = x2 - x1;
    long dy = y2 - y1;
    return std::sqrt(dx*dx + dy*dy);
}

void Util::splitAndCascade(std::vector<int>& first, std::vector<int>& second, int splitPoint)
{
    static thread_local std::vector<int> firstSplit;
    static thread_local std::vector<int> secondSplit;

    firstSplit.resize(first.size() - splitPoint);
    secondSplit.resize(second.size() - splitPoint);

    std::copy(first.begin() + splitPoint, first.end(), firstSplit.begin());
    std::copy(second.begin() + splitPoint, second.end(), secondSplit.begin());

    first.resize(splitPoint);
    second.resize(splitPoint);

    first.insert(first.end(), secondSplit.begin(), secondSplit.end());
    second.insert(second.end(), firstSplit.begin(), firstSplit.end());
}

void Util::splitAndFlipCascade(std::vector<int>& first, std::vector<int>& second, int splitPoint)
{
    static thread_local std::vector<int> firstSplit;
    static thread_local std::vector<int> secondSplit;

    firstSplit.resize(first.size() - splitPoint);
    secondSplit.resize(splitPoint);

    std::copy(first.begin() + splitPoint, first.end(), firstSplit.begin());
    std::copy(second.begin(), second.begin() + splitPoint, secondSplit.begin());

    first.resize(splitPoint);
    std::move(second.begin() + splitPoint, second.end(), second.begin());
    second.resize(second.size() - splitPoint);

    first.insert(first.end(), secondSplit.begin(), secondSplit.end());
    second.insert(second.end(), firstSplit.begin(), firstSplit.end());
}

}//cvrp namespace
