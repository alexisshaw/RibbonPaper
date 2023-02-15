//
// Created by 13114347 on 20/09/2020.
//

#include "util/getRandGen.h"

namespace SteaneSim{
    std::shared_ptr<std::mt19937_64> getRandGen() {
        thread_local std::shared_ptr<std::mt19937_64> rand_gen;
        thread_local bool is_initialised(false);
        if (!is_initialised) {
            std::random_device rd{"/dev/random"};
            std::seed_seq seq{rd(), rd(), rd(), rd()};
            rand_gen = std::make_shared<std::mt19937_64>(seq);
            is_initialised = true;
        }
        return rand_gen;
    }
}