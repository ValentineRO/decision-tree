#pragma once
#include "models.h"
#include "gurobi_c++.h"
using namespace std;

class solClust{
 public:
  vector<Tree> T;
  time_t totTime, clTime;
  int nbIter, nbClInitial, nbClFinal;
  double pseudoGap;
  vector<time_t> solveTime, centeringTime;
  vector<double> gap;
  vector<int> nbCl, obj, errCl, errDt;
  clustering finalCl;

  solClust(){
    T = {};
    totTime = 0;
    nbIter = 0;
    nbClInitial = 0;
    nbClFinal = 0;
    solveTime = {};
    centeringTime = {};
    gap = {};
    nbCl = {};
    obj = {};
    errCl = {};
    errDt = {};
  }

  void addClusteringTime(time_t clT){totTime += clT; clTime = clT;}
  void write(string filename);
};

bool isIntersecting(Tree T, clustering& cl, dataset& intialDt, dataset& currentDt, bool univ, int Nmin = 0);
solClust approxIteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, double timeL=3600);
solClust iteratingCART(dataset& dt, clustering& cl, parameters p);
solClust iteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, double timeL=3600);
/*
void buildPR(GRBModel& md, dataset& dt, clustering& cl, model_type modelt, parameters p);
solClust iteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, int timeL=3600);
*/
