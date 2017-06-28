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

void Util::swapElements(std::vector<int>& arr, int index1, int index2)
{
    int temp = arr[index1];
    arr[index1] = arr[index2];
    arr[index2] = temp;
}

void Util::randomShuffleElements(std::vector<int>& arr)
{
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(arr.begin(), arr.end(), g);
}

int Util::generateRandomNumberInRange(int min, int max)
{
    static thread_local std::default_random_engine gen;
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(gen);
}

void Util::splitAndCascade(std::vector<int>& first, std::vector<int>& second, int splitPoint)
{
    std::vector<int> firstSplit(first.begin() + splitPoint, first.end());
    std::vector<int> secondSplit(second.begin() + splitPoint, second.end());

    first.erase(first.begin() + splitPoint, first.end());
    second.erase(second.begin() + splitPoint, second.end());

    first.insert(first.end(), secondSplit.begin(), secondSplit.end());
    second.insert(second.end(), firstSplit.begin(), firstSplit.end());
}

void Util::splitAndFlipCascade(std::vector<int>& first, std::vector<int>& second, int splitPoint)
{
    std::vector<int> firstSplit(first.begin() + splitPoint, first.end());
    std::vector<int> secondSplit(second.begin(), second.begin() + splitPoint);

    first.erase(first.begin() + splitPoint, first.end());
    second.erase(second.begin(), second.begin() + splitPoint);

    first.insert(first.end(), secondSplit.begin(), secondSplit.end());
    second.insert(second.end(), firstSplit.begin(), firstSplit.end());
}

void Util::spliceAndCascade(std::vector<int>& first, std::vector<int>& second, int start, int end)
{
    std::vector<int> firstSplice(first.begin() + start, first.begin() + end);
    std::vector<int> secondSplice(second.begin() + start, second.begin() + end);

    first.erase(first.begin() + start, first.begin() + end);
    second.erase(second.begin() + start, second.begin() + end);

    first.insert(first.end(), secondSplice.begin(), secondSplice.end());
    second.insert(second.end(), firstSplice.begin(), firstSplice.end());
}

}//cvrp namespace
