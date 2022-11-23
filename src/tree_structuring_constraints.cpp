#include "tree_structuring_constraints.h"

int get_A(int t, int N){
  int A = (N + t + 1) / 2 - 1;
  if (t % 2 == 1){
    while ((A % 2 == 0)&&(A >= 1)){
      A = (A-1) / 2;
    }
    A = A / 2;
  }
  return A;
}


tree_structuring_constraints::tree_structuring_constraints(GRBModel& md, variables var, model_type mt, parameters p){
  string constraint_name;

  if (mt.base!=baseModel::F && mt.univ){
    for (int t=0; t<p.N; t++){ // sum a_jt = d_t
      GRBLinExpr sum_a_jt = 0;
      for (int j=0; j<p.J; j++){
	sum_a_jt += var.a[t*p.J+j];
      }
      constraint_name = "sumajt_EQ_dt_t=" + to_string(t);
      md.addConstr(sum_a_jt, GRB_EQUAL, var.d[t],constraint_name);
    }
    for (int t=0; t<p.N; t++){ // 0 <= b_t <= d_t
      constraint_name = "bt_LEQ_dt_t=" + to_string(t);
      md.addConstr(var.b[t], GRB_LESS_EQUAL, var.d[t], constraint_name);
      constraint_name = "bt_GEQ_0_t=" + to_string(t);
      md.addConstr(var.b[t], GRB_GREATER_EQUAL, 0, constraint_name);
    }
  }

  if (mt.base == baseModel::F){
     for (int t=0; t<p.N; t++){ // sum a_jt/d_t + sum c_kt = 1
        GRBLinExpr sum_a_jt = 0;
	if (mt.univ){
	  for (int j=0; j<p.J; j++){
            sum_a_jt += var.a[t*p.J+j];
	  }
	}
	else{
	  sum_a_jt += var.d[t];
	}
        GRBLinExpr sum_c_kt = 0;
        for (int k=0; k<p.K; k++){
            sum_c_kt += var.c[t*p.K+k];
        }
        constraint_name = "node_or_leaf_t=" + to_string(t);
        md.addConstr(sum_a_jt + sum_c_kt, GRB_EQUAL, 1, constraint_name);
    }
    for (int t=p.N; t<p.N+p.L; t++){ // sum c_kt = 1
      GRBLinExpr sum_c_kt = 0;
      for (int k=0; k<p.K; k++){
	sum_c_kt += var.c[t*p.K+k];
      }
      constraint_name = "one_classe_leaf_t=" + to_string(t);
      md.addConstr(sum_c_kt, GRB_EQUAL, 1, constraint_name);
    }
    if (mt.univ){
      for (int t=0; t<p.N; t++){ // b_t <= sum a_jt
	GRBLinExpr sum_a_jt = 0;
	for (int j=0; j<p.J; j++){
	  sum_a_jt += var.a[t*p.J+j];
	}
	constraint_name = "bt_domain_def_t="+to_string(t);
	md.addConstr(var.b[t],GRB_LESS_EQUAL,sum_a_jt,constraint_name);
      }
    }
  }
  
  if (mt.base!=baseModel::F || !mt.univ){
    for (int t = 1; t < p.N; t++) { // d_t <= d_a(t)
      int a_t = (t + 1) / 2 - 1;
      constraint_name = "dt_LEQ_dat_t=" + to_string(t);
      md.addConstr(var.d[t], GRB_LESS_EQUAL, var.d[a_t], constraint_name);
    }
  }
  
  if (!mt.univ){
    for (int t=0; t<p.N; t++){
      GRBLinExpr sum_s_jt = 0;
      for (int j=0; j<p.J; j++){
	constraint_name = "ajt_LEQ_ahjt_t=" + to_string(t) + "_j=" + to_string(j);
	md.addConstr(var.a[t*p.J+j], GRB_LESS_EQUAL, var.a_h[t*p.J+j], constraint_name);
	constraint_name = "-ajt_LEQ_ahjt_t=" + to_string(t) + "_j=" + to_string(j);
	md.addConstr(-var.a[t*p.J+j], GRB_LESS_EQUAL, var.a_h[t*p.J+j], constraint_name);
	constraint_name = "ajt_LEQ_sjt_t=" + to_string(t) + "_j=" + to_string(j);
	md.addConstr(var.a[t*p.J+j], GRB_LESS_EQUAL, var.s[t*p.J+j], constraint_name);
	constraint_name = "-ajt_LEQ_sjt_t=" + to_string(t) + "_j=" + to_string(j);
	md.addConstr(-var.a[t*p.J+j], GRB_LESS_EQUAL, var.s[t*p.J+j], constraint_name);
	constraint_name = "sjt_LEQ_dt=" + to_string(t) + "_j=" + to_string(j);
	md.addConstr(var.s[t*p.J+j], GRB_LESS_EQUAL, var.d[t], constraint_name);
	sum_s_jt += var.s[t*p.J+j];
      }
      constraint_name = "sum_sjt_GEQ_dt=" + to_string(t);
      md.addConstr(sum_s_jt,GRB_GREATER_EQUAL,var.d[t], constraint_name);
    }
    for (int t = 0; t < p.N; t++) { // sum ah_jt <= d_t
      GRBLinExpr sum_ah_jt = 0;
      for (int j = 0; j < p.J; j++) {
	sum_ah_jt += var.a_h[t * p.J + j];
      }
      constraint_name = "sumahjt_LEQ_dt_t=" + to_string(t);
      md.addConstr(sum_ah_jt, GRB_LESS_EQUAL, var.d[t], constraint_name);
    }
    /*
    for (int t = 0; t < p.N; t++) { // sum a_jt >= 0 (to avoid symetry)
      GRBLinExpr sum_a_jt = 0;
      for (int j = 0; j < p.J; j++) {
	sum_a_jt += var.a[t * p.J + j];
      }
      constraint_name = "sumajt_GEQ_0_t=" + to_string(t);
      md.addConstr(sum_a_jt, GRB_GREATER_EQUAL, 0, constraint_name);
    }
    */
    for (int t = 0; t < p.N; t++) { // -d_t <= b_t <= d_t
        constraint_name = "bt_LEQ_dt_t=" + to_string(t);
        md.addConstr(var.b[t], GRB_LESS_EQUAL, var.d[t], constraint_name);
	constraint_name = "bt_GEQ_-dt_t=" + to_string(t);
        md.addConstr(var.b[t], GRB_GREATER_EQUAL, -var.d[t], constraint_name);
    }
  }

  if (mt.eps){
    for (int t = 0; t < p.N; t++) {
      GRBLinExpr sum_ajt = 0;
      if (mt.base==baseModel::F && mt.univ){
	for (int j = 0; j < p.J; j++) {
	  sum_ajt += var.a[t * p.J + j];
	}
      }
      else{
	sum_ajt += var.d[t];
      }
      constraint_name = "epst_LEQ_sum_ajt_t=" + to_string(t);
      md.addConstr(var.eps[t], GRB_LESS_EQUAL, sum_ajt, constraint_name);
    }
  }
}
