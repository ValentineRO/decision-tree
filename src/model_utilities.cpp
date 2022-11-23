#include "model_utilities.h"

parameters::parameters(int depth, dataset &dt, double alp, int param_C, bool have_Lhat, int param_Nmin){
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

  alph = alp;
}

parameters parameters::parameters_copy(){
  parameters p = parameters();

  p.I = I;
  p.J = J;
  p.K = K;
  p.L_hat = L_hat;
  p.D = D;
  p.L = L;
  p.N = N;
  p.E = E;

  p.C = C;
  p.Nmin = Nmin;

  p.alph = alph;
  p.mu_min = mu_min;
  p.mu_max = mu_max;

  p.mu_vect = new double[J];
  for (int j=0; j<J; j++){
    p.mu_vect[j] = mu_vect[j];
  }

  return p;
}

void parameters::update(dataset &dt){
  I = dt.I;
  if (L_hat != 1){
    L_hat = dt.L_h;
  }
  mu_min = dt.mu_min;
  mu_max = dt.mu_max;

  for (int j=0; j<J; j++){
    mu_vect[j] = dt.mu_vect[j];
  }
  
}

variable_def::variable_def() {
    lb = {};
    ub = {};
    coef = {};
    type = {};
    names = {};
    count = 0;
}

variable_def::variable_def(string name, int nb_var, bool isBinary, double LB, double UB) {
  count = nb_var;
  lb = new double[nb_var];
  ub = new double[nb_var];
  coef = new double[nb_var];
  type = new char[nb_var];
  names = new string[nb_var];
  for (int i = 0; i < nb_var; i++) {
    lb[i] = LB;
    ub[i] = UB;
    coef[i] = 0.0;
    if (isBinary){
      type[i] = GRB_BINARY;
    }
    else{
      type[i] = GRB_CONTINUOUS;
    }
    names[i] = name + to_string(i);
  }
}

variables::variables(GRBModel& md, model_type mt, parameters p){
  variable_def a_var;
  if (mt.univ){
    a_var = variable_def("a", p.J * p.N, !mt.relaxation);
  }
  else{
    a_var = variable_def("a", p.J * p.N, false, -1, 1);
  }
  a = md.addVars(a_var.lb, a_var.ub, a_var.coef, a_var.type, a_var.names, a_var.count);

  if (!mt.univ){
    variable_def ah_var = variable_def("ah", p.J * p.N, false),
      s_var = variable_def("s", p.J * p.N, true);
    s = md.addVars(s_var.lb, s_var.ub, s_var.coef, s_var.type, s_var.names, s_var.count);
    a_h = md.addVars(ah_var.lb, ah_var.ub, ah_var.coef, ah_var.type, ah_var.names, ah_var.count);
  }

  variable_def b_var;
  if (mt.univ){
    b_var = variable_def("b", p.N, false);    
  }
  else{
    b_var = variable_def("b", p.N, false, -1, 1);
  }
  b = md.addVars(b_var.lb, b_var.ub, b_var.coef, b_var.type, b_var.names, b_var.count);

  variable_def c_var;
  if (mt.base == baseModel::F){
    c_var = variable_def("c", p.K * (p.L+p.N), !mt.relaxation);
  }
  else{
    c_var = variable_def("c", p.K * p.L, !mt.relaxation);
  }
  c = md.addVars(c_var.lb, c_var.ub, c_var.coef, c_var.type, c_var.names, c_var.count);

  if ((!mt.univ) or (mt.base != baseModel::F)){
    variable_def d_var = variable_def("d", p.N, !mt.relaxation);
    d = md.addVars(d_var.lb, d_var.ub, d_var.coef, d_var.type, d_var.names, d_var.count);
  }

  if (mt.eps){
    variable_def eps_var = variable_def("eps", p.N, false);
    eps = md.addVars(eps_var.lb, eps_var.ub, eps_var.coef, eps_var.type, eps_var.names, eps_var.count);
  }

  if (mt.base != baseModel::F){
    variable_def l_var = variable_def("l", p.L, !mt.relaxation);
    l = md.addVars(l_var.lb, l_var.ub, l_var.coef, l_var.type, l_var.names, l_var.count);
  }

  if (mt.base == baseModel::OCT){
    variable_def Lt_var = variable_def("Lt", p.L, false, 0, p.I),
      Nt_var = variable_def("Nt", p.L, false, 0, p.I),
      Nkt_var = variable_def("Nkt", p.L * p.K, false, 0, p.I);
    Lt = md.addVars(Lt_var.lb, Lt_var.ub, Lt_var.coef, Lt_var.type, Lt_var.names, Lt_var.count);
    Nt = md.addVars(Nt_var.lb, Nt_var.ub, Nt_var.coef, Nt_var.type, Nt_var.names, Nt_var.count);
    Nkt = md.addVars(Nkt_var.lb, Nkt_var.ub, Nkt_var.coef, Nkt_var.type, Nkt_var.names, Nkt_var.count);
  }

  if (mt.base == baseModel::FOCT){
    variable_def theta_var = variable_def("theta", p.L * p.I * p.K, false);
    theta = md.addVars(theta_var.lb, theta_var.ub, theta_var.coef, theta_var.type, theta_var.names, theta_var.count);
  }
  if (mt.base == baseModel::GOCT){
    variable_def theta_var = variable_def("theta", p.L * p.K, false, 0, p.I);
    theta = md.addVars(theta_var.lb, theta_var.ub, theta_var.coef, theta_var.type, theta_var.names, theta_var.count);
  }

  if (mt.base == baseModel::F){
    variable_def u_var = variable_def("u", p.I * p.E, !mt.relaxation);
    u = md.addVars(u_var.lb, u_var.ub, u_var.coef, u_var.type, u_var.names, u_var.count);
  }
  else{
    variable_def z_var = variable_def("z", p.I * p.L, !mt.relaxation);
    z = md.addVars(z_var.lb, z_var.ub, z_var.coef, z_var.type, z_var.names, z_var.count);
  }
}
