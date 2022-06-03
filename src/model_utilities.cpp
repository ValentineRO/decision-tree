#include "model_utilities.h"

binary_variable_def::binary_variable_def() {
    lb = {};
    ub = {};
    coef = {};
    type = {};
    names = {};
    count = 0;
}

binary_variable_def::binary_variable_def(string name, int nb_var, double UB) {
  count = nb_var;
  lb = new double[nb_var];
  ub = new double[nb_var];
  coef = new double[nb_var];
  type = new char[nb_var];
  names = new string[nb_var];
  for (int i = 0; i < nb_var; i++) {
    lb[i] = 0.0;
    ub[i] = UB;
    coef[i] = 0.0;
    type[i] = GRB_CONTINUOUS;
    names[i] = name + to_string(i);
  }
}

parameters::parameters(int depth, dataset &dt, double alp, double param_mu, int param_C, bool have_Lhat, int param_Nmin){
  I = dt.I;
  J = dt.J;
  K = dt.K;
  
  if (have_Lhat){
    L_hat = dt.L_h;
  }
  else{
    L_hat = 1;
  }

  mu_min = dt.mu_min;
  mu_max = dt.mu_max;
  mu_vect = dt.mu_vect;

  D = depth;
  L = pow(2,D);
  N = pow(2,D) - 1;
  E = 3*N + L + 1;

  Nmin = param_Nmin;
  C = param_C;

  mu = param_mu;
  alph = alp;
}

variables::variables(GRBModel& md, model_type mt, parameters p){
  if (mt.relaxation){
    if (mt.univ){
      binary_variable_def a_var = binary_variable_def("a", p.J * p.N);
      a = md.addVars(a_var.lb, a_var.ub, a_var.coef, a_var.type, a_var.names, a_var.count);
    }
    else{
      binary_variable_def s_var = binary_variable_def("s", p.J * p.N);
      s = md.addVars(s_var.lb, s_var.ub, s_var.coef, s_var.type, s_var.names, s_var.count);
    }
    if (mt.base == baseModel::F){
      binary_variable_def c_var = binary_variable_def("c", p.K * (p.L+p.N)),
	u_var = binary_variable_def("u", p.E*p.I);
      
      c = md.addVars(c_var.lb, c_var.ub, c_var.coef, c_var.type, c_var.names, c_var.count);
      u = md.addVars(u_var.lb, u_var.ub, u_var.coef, u_var.type, u_var.names, u_var.count);
      if (!mt.univ){
	binary_variable_def d_var = binary_variable_def("d", p.N);
	d = md.addVars(d_var.lb, d_var.ub, d_var.coef, d_var.type, d_var.names, d_var.count);
      }
    }
    else{
      binary_variable_def c_var = binary_variable_def("c", p.K * (p.L+p.N)),
	d_var = binary_variable_def("d", p.N),
	l_var = binary_variable_def("l", p.L),
	z_var = binary_variable_def("z", p.I*p.L);
      c = md.addVars(c_var.lb, c_var.ub, c_var.coef, c_var.type, c_var.names, c_var.count);
      d = md.addVars(d_var.lb, d_var.ub, d_var.coef, d_var.type, d_var.names, d_var.count);
      l = md.addVars(l_var.lb, l_var.ub, l_var.coef, l_var.type, l_var.names, l_var.count);
      z = md.addVars(z_var.lb, z_var.ub, z_var.coef, z_var.type, z_var.names, z_var.count);
    }
  }
  else{
    if (mt.univ){
      a = md.addVars(p.J * p.N, GRB_BINARY);
    }
    else{
      s = md.addVars(p.J * p.N, GRB_BINARY);
    }
    if (mt.base == baseModel::F){
      c = md.addVars((p.N+p.L)*p.K,GRB_BINARY);
      u = md.addVars(p.E*p.I, GRB_BINARY);
      if (!mt.univ){
	d = md.addVars(p.N, GRB_BINARY);
      }
    }
    else{
      c = md.addVars(p.L*p.K, GRB_BINARY);
      d = md.addVars(p.N, GRB_BINARY);
      l = md.addVars(p.L, GRB_BINARY);
      z = md.addVars(p.I * p.L, GRB_BINARY);
    }
  }
  
  b = md.addVars(p.N, GRB_CONTINUOUS);
  if (!mt.univ){
    a = md.addVars(p.J * p.N, GRB_CONTINUOUS);
    a_h = md.addVars(p.J * p.N, GRB_CONTINUOUS);
  }
  if (mt.eps){
    eps = md.addVars(p.N, GRB_CONTINUOUS);
  }

  if (mt.base == baseModel::OCT){
    Lt = md.addVars(p.L, GRB_CONTINUOUS);
    Nt = md.addVars(p.L, GRB_CONTINUOUS);
    Nkt = md.addVars(p.L*p.K, GRB_CONTINUOUS);
  }

  if (mt.base == baseModel::FOCT){
    theta = md.addVars(p.L*p.I*p.K, GRB_CONTINUOUS);
  }

  if (mt.base == baseModel::GOCT){
    theta = md.addVars(p.L*p.K, GRB_CONTINUOUS);
  }
  
}
