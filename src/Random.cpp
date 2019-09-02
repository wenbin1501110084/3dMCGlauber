// Copyright @ Chun Shen 2018

#include "Random.h"

namespace RandomUtil {

Random::Random(int seed_in, double min, double max, int seed_add) :
    rand_uniform_dist(min, max) {
    seed = seed_in;
    if (seed == -1) {
        seed = ran_dev() + seed_add;
    }
    ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seed));
}

}

