#pragma once
#include "model_utilities.h"
#include "gurobi_c++.h"
#include <string>
#include <vector>
using namespace std;

class following_constraints{
 public:
  following_constraints() {}
  following_constraints(GRBModel& md, variables var, model_type mt, parameters p, dataset& dt);
};
