#pragma once
#include "iteratingOTP.h"
#include "gurobi_c++.h"
using namespace std;

class training_results{
 public:
  double globalTime,
    globalNbOpti,
    globalNbIter;
  // if there is post-processing we will have 3 differents outcomes : with, without and both (best validation score is selected)
  double* time;
  double* nbOpti;
  double* nbIter;
  double* percErrWithoutPP;
  double* percErrWithPP;
  double* percErrBestOfBoth;
  double* alphWithoutPP;
  double* alphWithPP;
  double* alphBestOfBoth;
  string* bestTreeWithoutPP;
  string* bestTreeWithPP;
  string* bestTreeBestOfBoth;

  training_results(){}

  training_results(int nbDiffDepth){
    globalTime = 0;
    globalNbOpti = 0;
    globalNbIter = 0;
    time = new double[nbDiffDepth];
    nbOpti = new double[nbDiffDepth];
    nbIter = new double[nbDiffDepth];
    percErrWithoutPP = new double[nbDiffDepth];
    percErrWithPP = new double[nbDiffDepth];
    percErrBestOfBoth = new double[nbDiffDepth];
    alphWithoutPP = new double[nbDiffDepth];
    alphWithPP = new double[nbDiffDepth];
    alphBestOfBoth = new double[nbDiffDepth];
    bestTreeWithoutPP = new string[nbDiffDepth]; 
    bestTreeWithPP = new string[nbDiffDepth]; 
    bestTreeBestOfBoth = new string[nbDiffDepth]; 
  }
};

class training_results_cl{
 public :
  double globalTime,
    globalNbOpti,
    globalNbIter,
    globalNbIterCl,
    globalPseudoGap,
    globalFinalReduction;

  double* time;
  double* nbIter;
  double* nbOpti;
  double* nbIterCl;
  double* pseudoGap;
  double* finalReduction;
  double* errTrain;
  double* errTest;
  string* bestTree;
  
  training_results_cl(){}

  training_results_cl(int nbDiffDepth){
    globalTime = 0;
    globalNbOpti = 0;
    globalNbIter = 0;
    globalNbIterCl = 0;
    globalPseudoGap = 0;
    globalFinalReduction = 0;
    time = new double[nbDiffDepth];
    nbIter = new double[nbDiffDepth];
    nbOpti = new double[nbDiffDepth];
    nbIterCl = new double[nbDiffDepth];
    pseudoGap = new double[nbDiffDepth];
    finalReduction = new double[nbDiffDepth];
    errTrain = new double[nbDiffDepth];
    errTest = new double[nbDiffDepth];
    bestTree = new string[nbDiffDepth];
  }
};

training_results learning_Bertsimas(dataset& dt_train, dataset& dt_validation, dataset& dt_test, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin);

training_results learning(dataset& dt_train, dataset& dt_validation, dataset& dt_test, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin);

training_results_cl learningWithClustering(dataset& dt_train, dataset& dt_validation, dataset& dt_test, clustering cl, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin, int selectingStrat);
