#pragma once
#include <string>
#include "solution.h"
#include "gurobi_c++.h"
using namespace std;

enum baseModel {
  F,
  OCT,
  QOCT,
  FOCT,
  GOCT
};

class model_type{
 public :
  baseModel base; // F, OCT, QOCT, FOCT, GOCT
  bool univ;
  bool eps;
  bool C_const;
  bool relaxation;

  model_type(){};
  model_type(baseModel bm, bool u, bool e, bool c, bool r){
    base = bm;
    univ = u;
    eps = e;
    C_const = c;
    relaxation = r;
  }
};

class parameters{
 public :
  int I,J,K,L_hat;
  int D,L,N,E;

  int Nmin,C;

  double alph,mu_min,mu_max;
  double* mu_vect;

  parameters(){};
  parameters(int depth, dataset &dt, double alp, int param_C = -1, bool have_Lhat=true, int param_Nmin = 0);

  parameters parameters_copy();
  void update(dataset &dt);
};


class variable_def {
public:
    double* lb;
    double* ub;
    double* coef;
    char* type;
    string* names;
    int count;

    variable_def();
    variable_def(string name, int nb_var, bool isBinary, double LB=0.0, double UB=1.0);
};

class variables {
 public :
  GRBVar* a; // a_jt = a[t*J + j]
  GRBVar* a_h;
  GRBVar* s;
  GRBVar* b;
  GRBVar* c; // c_kt = c[t*K + k]
  GRBVar* d;
  GRBVar* eps;
  GRBVar* l;
  GRBVar* Lt;
  GRBVar* Nt;
  GRBVar* Nkt;
  GRBVar* theta; // Theta_kt = Theta[t*K + k] et theta_tik = theta[t*I*K + i*K + k]
  GRBVar* u; // u_ie = [i*E + e]
  GRBVar* z; // z_it = [i*N + t]

  variables() {};
  variables(GRBModel& md, model_type mt, parameters p);
};
