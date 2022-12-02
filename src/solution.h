#pragma once
#include "tree.h"
#include <utility>
// #include <vector>
using namespace std;

class solution{
public :
    Tree T;
    double obj;
    vector<double> objEvo;
    int error_train;
    int nb_br;
    float time;
    double gap;
    int nodes;
    double root_rel;

    solution() {}
    solution(Tree T_, double o, vector<double> oEvo, int et, int nb, float t, double g, int n, double r);
};

class optimal_tree{
 public:
  int D,C,E,Ev,Et;
  string namefile;
  
  optimal_tree(){}
  optimal_tree(int d, int c, int e, int ev, int et, string nf){
    D = d;
    C = c;
    E = e;
    Ev = ev;
    Et = et;
    namefile = nf;
  }
};

class dominating_trees{
 public:
  optimal_tree* trees;
  int Lh;
  double* alph;
  int nb_sol;

  dominating_trees(){};
  dominating_trees(int nb_sol_tot, int lh){
    trees = new optimal_tree[nb_sol_tot];
    alph = new double[nb_sol_tot+1];
    nb_sol = 0;
    Lh = lh;
  }

  bool add_tree(optimal_tree& opt);

  int BestTree(bool smallest_alpha);
  int BestWarmstart(int D, int C);
};
