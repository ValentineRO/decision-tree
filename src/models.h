#pragma once
#include "objective.h"
#include "branching_constraints.h"
#include "following_constraints.h"
#include "tree_structuring_constraints.h"
#include "counting_errors.h"
#include "clustering.h"
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;

class build_model{
 public :
  model_type mt;
  parameters param;
  variables var;

  objective obj;
  following_constraints flwc;
  tree_structuring_constraints tsc;
  branching_constraints bc;
  counting_errors ce;

  build_model() {}
  build_model(GRBModel& md, dataset& dt, model_type modelt, parameters param);

  void add_warmstart(Tree T, dataset& dt);
  void addWarmStartInConstraints(GRBModel& md, dataset& dt, Tree T);
  // void addClusteringWarmStart();

  double read_root_rel_objvalue();
  int compute_number_branchings();
  void buildTree(Tree& T);
  void get_z(int z[]);
  
  solution solve(GRBModel& md, double timeL = 3600.0, int solLim = -1);
  solution solve2(GRBModel& md, double timeL = 3600.0, double stoppingPerc = 0.333);
};

