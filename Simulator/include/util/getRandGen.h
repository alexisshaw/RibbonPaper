//
// Created by 13114347 on 20/09/2020.
//

#ifndef STEANE_SIM_MYRANDGEN_H
#define STEANE_SIM_MYRANDGEN_H
#include <random>
#include <memory>

namespace SteaneSim{
    std::shared_ptr<std::mt19937_64> getRandGen();
}


#endif //STEANE_SIM_MYRANDGEN_H
