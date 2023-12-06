#pragma once
#include "gurobi_c++.h"

class GurobiEnvironment {
public:
    static GurobiEnvironment& getInstance();
    GRBEnv& getEnvironment();

private:
    GurobiEnvironment(); // Private constructor

    GRBEnv env;
};

extern GurobiEnvironment& gurobiEnv;
