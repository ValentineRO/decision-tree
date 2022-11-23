#pragma once
#include "models.h"
#include "gurobi_c++.h"
using namespace std;

class solClust{
 public:
  Tree T;
  time_t totTime, clTime;
  int nbIter, nbClInitial, nbClFinal, errTr, errTst;
  vector<time_t> solveTime, centeringTime;
  vector<bool> isOpti;
  vector<int> nbCl, errCl, errDt;
  clustering finalCl;

  solClust(){
    T = Tree();
    totTime = 0;
    nbIter = 0;
    errTr = 0;
    errTst = 0;
    nbClInitial = 0;
    nbClFinal = 0;
    solveTime = {};
    centeringTime = {};
    isOpti = {};
    nbCl = {};
    errCl = {};
    errDt = {};
  }

  void addClusteringTime(time_t clT){totTime += clT; clTime = clT;}
  void write(string filename);
};

bool isIntersecting(Tree T, clustering& cl, dataset& intialDt, dataset& currentDt, bool univ, bool updateCl);
solClust iteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, int timeL=3600);
