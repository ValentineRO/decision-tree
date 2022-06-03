#pragma once
#include "model_utilities.h"
#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
using namespace std;

class objective{
 public :
  objective(){}
  objective(GRBModel& md, variables var, model_type mt, parameters p, dataset& dt);
};
