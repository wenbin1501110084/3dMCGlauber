// Copyright (C) 2018 Chun Shen
#include "ParametersMap.h"
#include "Util.h"
#include <fstream>
#include <iostream>

ParametersMap::~ParametersMap() {
    parameter_map.clear();
}


void ParametersMap::read_in_parameters_from_file(string filename) {
    std::ifstream input_file(filename);
    if (!input_file.is_open()) {
        std::cout << "[ParametersMap]: can not open the file: " << filename 
                  << std::endl;
        exit(1);
    }
    string line;
    while (getline(input_file,line)) {
        auto param = StringUtility::parse_a_line(line, " ", "#");
        set_parameter(param[0], param[1]);
    }
}
