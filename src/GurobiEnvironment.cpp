#include "GurobiEnvironment.h"

// Define the methods

GurobiEnvironment::GurobiEnvironment() {
    // Initialize the Gurobi environment with desired settings
    env.set(GRB_IntParam_Threads, 2);
}

GurobiEnvironment& GurobiEnvironment::getInstance() {
    static GurobiEnvironment instance;
    return instance;
}

GRBEnv& GurobiEnvironment::getEnvironment() {
    return env;
}

GurobiEnvironment& gurobiEnv = GurobiEnvironment::getInstance();
