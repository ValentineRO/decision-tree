#pragma once
#include "gurobi_c++.h"
#include <string>
#include "model_utilities.h"
using namespace std;

class counting_errors{
public :
    counting_errors() {};
    counting_errors(GRBModel& md, variables var, model_type mt, parameters p, dataset& dt);
};
