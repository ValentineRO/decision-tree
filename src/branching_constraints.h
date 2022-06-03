#pragma once
#include <math.h>
#include <vector>
#include <string>
#include "gurobi_c++.h"
#include "model_utilities.h"
using namespace std;

void get_ancestors(vector<int>& A_L, vector<int>& A_R, int t, int D);

class branching_constraints{
 public:
  branching_constraints() {};
  branching_constraints(GRBModel& md, variables var, model_type mt, parameters p, dataset& dt);
};
