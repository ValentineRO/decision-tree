#pragma once
#include "models.h"
#include "gurobi_c++.h"
using namespace std;

class solClust{
 public:
  Tree T;
  time_t totTime;
  int nbIter, nbClInitial, nbClFinal, errTr, errTst;
  vector<time_t> solveTime, centeringTime;
  vector<bool> isOpti;
  vector<int> nbCl, errCl, errDt;

  solClust(){
    T = Tree();
    totTime = 0;
    nbIter, errTr, errTst, nbClInitial, nbClFinal = 0, 0, 0, 0, 0;
    solveTime = {};
    centeringTime = {};
    isOpti = {};
    nbCl = {};
    errCl = {};
    errDt = {};
  }

  void write(string filename);
};

bool isIntersecting(Tree T, clustering& cl, dataset& intialDt, dataset& currentDt, bool univ, bool updateCl);
solClust decontractingAlgorithm(dataset& dt, clustering& cl, model_type modelt, parameters p, int timeL=3600);
