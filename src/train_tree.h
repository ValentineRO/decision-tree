#pragma once
#include "iteratingOTP.h"
#include "gurobi_c++.h"
using namespace std;

class training_results{
 public:
  double time,
    nb_opti,
    nb_iter;
  // if there is post-processing we will have 3 differents outcomes : with, without and both (best validation score is selected)
  double* perc_err;
  double* alph;
  string* best_tree;

  training_results(){
    time = 0;
    nb_opti = 0;
    nb_iter = 0;
    perc_err = new double[6];
    alph = new double[6];
    best_tree = new string[6];
    for (int i=0; i<6; i++){
      perc_err[i] = 0;
      alph[i] = 0;
    }
  }
};

training_results learning_Bertsimas(dataset& dt_train, dataset& dt_validation, dataset& dt_test, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin);

training_results learning(dataset& dt_train, dataset& dt_validation, dataset& dt_test, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin);
