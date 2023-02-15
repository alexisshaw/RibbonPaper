
#include <random>
#include <iostream>
#include "sim/norm.h"

unsigned int random_seed(){
    //std::cerr << "init rng" << std::endl;
    return std::random_device{}();
}

double ds_uniform()
{
    thread_local std::mt19937_64 generator(random_seed());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
//return 0.5;
}

double ds_norm()
{
    thread_local std::mt19937_64 generator(random_seed());
    std::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
//return 0.0;
}
