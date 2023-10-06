#include "branching_constraints.h"

void get_ancestors(vector<int>& A_L, vector<int>& A_R, int t, int D){
    int t1 = t,
        t2 = t + (int)pow(2,D);
    for (int d = 0; d < D; d++) {
        t2 = t2 / 2;
        if (t1 % 2 == 0){
            A_L.push_back(t2-1);
        }
        else{
            A_R.push_back(t2-1);
        }
        t1 = t1 / 2;
    }
}

branching_constraints::branching_constraints(GRBModel& md, variables var, model_type mt, parameters p, dataset& dt){
  string constraint_name;
  if (mt.base == baseModel::F){
    for (int t=0; t<p.N; t++){
      for (int i=0; i<p.I; i++){
	int index_u_tlt = i * p.E + t * 3 + 1,
	  index_u_trt = i * p.E + t * 3 + 2;
	
	GRBLinExpr left_b_mb_g = 0,
	  left_b_mb_d = 0;
	for (int j = 0; j < p.J; j++) {
	  left_b_mb_g += var.a[t * p.J + j] * dt.X[i * p.J + j];
	  if (mt.univ){
	    left_b_mb_g += var.a[t * p.J + j] * (p.mu_vect[j] - p.mu_min);
	  }
	}
	if (mt.univ){
	  left_b_mb_g += p.mu_min;
	  left_b_mb_d += (1+p.mu_max)*(1-var.u[index_u_tlt]);
	}
	else{
	  left_b_mb_g += mu;
	  left_b_mb_d += (2+mu)*(1-var.u[index_u_tlt]);
	}
	if (mt.eps){
	  left_b_mb_g += var.eps[t];
	}
	left_b_mb_d += var.b[t];

	constraint_name = "left_branching_t=" + to_string(t) + "_i=" + to_string(i);
	md.addConstr(left_b_mb_g, GRB_LESS_EQUAL, left_b_mb_d, constraint_name);
	
	GRBLinExpr right_b_mb_g = 0,
	    right_b_mb_d = 0;
	for (int j = 0; j < p.J; j++) {
	  right_b_mb_g += var.a[t * p.J + j] * dt.X[i * p.J + j];
	}
	if (mt.univ){
	  right_b_mb_d -= 1-var.u[index_u_trt];
	}
	else{
	  right_b_mb_d -= 2*(1-var.u[index_u_trt]);
	}
	right_b_mb_d += var.b[t];
	constraint_name = "right_branching_t=" + to_string(t) + "_i=" + to_string(i);
	md.addConstr(right_b_mb_g, GRB_GREATER_EQUAL, right_b_mb_d, constraint_name);
      }
    }
  }
  else{
    for (int t=0; t<p.L; t++){
      vector<int> A_L = {};
      vector<int> A_R = {};
      get_ancestors(A_L,A_R,t,p.D);
      for (int i=0; i<p.I; i++){
	for (auto g : A_L) {
	  GRBLinExpr left_b_mb_g = 0,
	    left_b_mb_d = 0;
	  for (int j = 0; j < p.J; j++) {
	    left_b_mb_g += var.a[g * p.J + j] * dt.X[i * p.J + j];
	    if (mt.univ){
	      left_b_mb_g += var.a[g * p.J + j] * (p.mu_vect[j] - p.mu_min);
	    }
	  }
	  if (mt.univ){
	    left_b_mb_g += p.mu_min;
	    left_b_mb_d += (1+p.mu_max)*(1-var.z[i*p.L+t]);
	  }
	  else{
	    left_b_mb_g += mu;
	    left_b_mb_d += (2+mu)*(1-var.z[i*p.L+t]);
	  }
	  if (mt.eps){
	    left_b_mb_g += var.eps[g];
	  }
	  left_b_mb_d += var.b[g];
                
	  constraint_name = "left_branching_t=" + to_string(t) + "_i=" + to_string(i) + "_m=" + to_string(g);
	  md.addConstr(left_b_mb_g,GRB_LESS_EQUAL,left_b_mb_d,constraint_name);
	}
	for (auto g : A_R){
	  GRBLinExpr right_b_mb_g = 0,
	    right_b_mb_d = 0;
	  for (int j = 0; j < p.J; j++) {
	    right_b_mb_g += var.a[g * p.J + j] * dt.X[i * p.J + j];
	  }
	  if (mt.univ){
	    right_b_mb_d -= 1-var.z[i*p.L+t];
	  }
	  else{
	    right_b_mb_d -= 2*(1-var.z[i*p.L+t]);
	  }
	  right_b_mb_d += var.b[g];
	  constraint_name = "right_branching_t=" + to_string(t) + "_i=" + to_string(i) + "_m=" + to_string(g);
	  md.addConstr(right_b_mb_g,GRB_GREATER_EQUAL,right_b_mb_d, constraint_name);
	}
      }
      if (mt.base!=baseModel::OCT){
	for (auto g: A_L){
	  constraint_name = "lt_LEQ_dm_t=" + to_string(t) + "_m=" + to_string(g);
	  md.addConstr(var.l[t],GRB_LESS_EQUAL,var.d[g], constraint_name);
	}
      }
    }
  }
}

branching_constraints::branching_constraints(GRBModel& md, variables var, model_type mt, parameters p, clustering& cl){
  string constraint_name;
  if (mt.base == baseModel::F){
    for (int t=0; t<p.N; t++){
      for (int i=0; i<cl.refDt->I; i++){
	int index_u_tlt = cl.placeOf[cl.clusterOf[i]] * p.E + t * 3 + 1,
	  index_u_trt = cl.placeOf[cl.clusterOf[i]] * p.E + t * 3 + 2;
	
	GRBLinExpr left_b_mb_g = 0,
	  left_b_mb_d = 0;
	for (int j = 0; j < p.J; j++) {
	  left_b_mb_g += var.a[t * p.J + j] * cl.refDt->X[i * p.J + j];
	  if (mt.univ){
	    left_b_mb_g += var.a[t * p.J + j] * (p.mu_vect[j] - p.mu_min);
	  }
	}
	if (mt.univ){
	  left_b_mb_g += p.mu_min;
	  left_b_mb_d += (1+p.mu_max)*(2-var.u[index_u_tlt]-var.r[i]);
	}
	else{
	  left_b_mb_g += mu;
	  left_b_mb_d += (2+mu)*(2-var.u[index_u_tlt]-var.r[i]);
	}
	if (mt.eps){
	  left_b_mb_g += var.eps[t];
	}
	left_b_mb_d += var.b[t];

	constraint_name = "left_branching_t=" + to_string(t) + "_i=" + to_string(i) + "_g=" + to_string(cl.placeOf[cl.clusterOf[i]]);
	md.addConstr(left_b_mb_g, GRB_LESS_EQUAL, left_b_mb_d, constraint_name);
	
	GRBLinExpr right_b_mb_g = 0,
	    right_b_mb_d = 0;
	for (int j = 0; j < p.J; j++) {
	  right_b_mb_g += var.a[t * p.J + j] * cl.refDt->X[i * p.J + j];
	}
	if (mt.univ){
	  right_b_mb_d -= 2-var.u[index_u_trt]-var.r[i];
	}
	else{
	  right_b_mb_d -= 2*(2-var.u[index_u_trt]-var.r[i]);
	}
	right_b_mb_d += var.b[t];
	constraint_name = "right_branching_t=" + to_string(t) + "_i=" + to_string(i) + "_g=" + to_string(cl.placeOf[cl.clusterOf[i]]);
	md.addConstr(right_b_mb_g, GRB_GREATER_EQUAL, right_b_mb_d, constraint_name);
      }
    }
  }
  else{
    for (int t=0; t<p.L; t++){
      vector<int> A_L = {};
      vector<int> A_R = {};
      get_ancestors(A_L,A_R,t,p.D);
      for (int i=0; i<cl.refDt->I; i++){
	for (auto g : A_L) {
	  GRBLinExpr left_b_mb_g = 0,
	    left_b_mb_d = 0;
	  for (int j = 0; j < p.J; j++) {
	    left_b_mb_g += var.a[g * p.J + j] * cl.refDt->X[i * p.J + j];
	    if (mt.univ){
	      left_b_mb_g += var.a[g * p.J + j] * (p.mu_vect[j] - p.mu_min);
	    }
	  }
	  if (mt.univ){
	    left_b_mb_g += p.mu_min;
	    left_b_mb_d += (1+p.mu_max)*(2-var.z[cl.placeOf[cl.clusterOf[i]]*p.L+t]-var.r[i]);
	  }
	  else{
	    left_b_mb_g += mu;
	    left_b_mb_d += (2+mu)*(2-var.z[cl.placeOf[cl.clusterOf[i]]*p.L+t]-var.r[i]);
	  }
	  if (mt.eps){
	    left_b_mb_g += var.eps[g];
	  }
	  left_b_mb_d += var.b[g];
                
	  constraint_name = "left_branching_t=" + to_string(t) + "_i=" + to_string(i) + "_g=" + to_string(cl.placeOf[cl.clusterOf[i]]) + "_m=" + to_string(g);
	  md.addConstr(left_b_mb_g,GRB_LESS_EQUAL,left_b_mb_d,constraint_name);
	}
	for (auto g : A_R){
	  GRBLinExpr right_b_mb_g = 0,
	    right_b_mb_d = 0;
	  for (int j = 0; j < p.J; j++) {
	    right_b_mb_g += var.a[g * p.J + j] * cl.refDt->X[i * p.J + j];
	  }
	  if (mt.univ){
	    right_b_mb_d -= 2-var.z[cl.placeOf[cl.clusterOf[i]]*p.L+t]-var.r[i];
	  }
	  else{
	    right_b_mb_d -= 2*(2-var.z[cl.placeOf[cl.clusterOf[i]]*p.L+t]-var.r[i]);
	  }
	  right_b_mb_d += var.b[g];
	  constraint_name = "right_branching_t=" + to_string(t) + "_i=" + to_string(i) + "_g=" + to_string(cl.placeOf[cl.clusterOf[i]]) + "_m=" + to_string(g);
	  md.addConstr(right_b_mb_g,GRB_GREATER_EQUAL,right_b_mb_d, constraint_name);
	}
      }
      if (mt.base!=baseModel::OCT){
	for (auto g: A_L){
	  constraint_name = "lt_LEQ_dm_t=" + to_string(t) + "_m=" + to_string(g);
	  md.addConstr(var.l[t],GRB_LESS_EQUAL,var.d[g], constraint_name);
	}
      }
    }
  }
  for (int g=0; g<cl.clusters.size(); g++){
    GRBLinExpr sum_r = 0;
    for (auto i: cl.clusters[g].pts){
      sum_r += var.r[i];
    }
    constraint_name = "sum_ri_=_1_g=" + to_string(g);
    md.addConstr(sum_r,GRB_EQUAL,1, constraint_name);
  }
}
