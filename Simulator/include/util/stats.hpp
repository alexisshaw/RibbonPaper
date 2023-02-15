//
// Created by 13114347 on 27/05/2020.
//

#ifndef STEANE_SIM_STATS_HPP
#define STEANE_SIM_STATS_HPP
#include <vector>
#include <numeric>

template<typename T>
double variance(const std::vector<T> &vec)
{
    size_t sz = vec.size();
    if (sz == 1)
        return 0.0;

    // Calculate the mean
    double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

    // Now calculate the variance
    auto variance_func = [mean, sz](T accumulator, const T& val)
    {
        return accumulator + ((val - mean)*(val - mean) / (sz - 1));
    };

    return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

template<typename T>
double standarddeviation(std::vector<T> samples)
{
    return sqrt(variance(samples));
}

template<typename T>
T median(std::vector<T> &v)
{
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}

template<typename T>
double mean(const std::vector<T> &vec)
{
    size_t sz = vec.size();
    if (sz == 1)
        return 0.0;

    // Calculate the mean
    double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

    return mean;
}
#endif //STEANE_SIM_STATS_HPP