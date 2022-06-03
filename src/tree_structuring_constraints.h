#pragma once
#include "gurobi_c++.h"
#include <string>
#include "model_utilities.h"
using namespace std;

int get_A(int t, int N);

class tree_structuring_constraints{
 public :
  tree_structuring_constraints() {}
  tree_structuring_constraints(GRBModel& md, variables var, model_type mt, parameters p);
};
