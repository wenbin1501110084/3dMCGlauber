// Copyright (C) 2018 Chun Shen

#include <sstream>
#include "Util.h"
#include "Parameters.h"

namespace MCGlb {

void Parameters::set_b(real b_in) {
    std::ostringstream param_val;
    param_val << b_in;
    set_parameter("b", param_val.str());
}

}

