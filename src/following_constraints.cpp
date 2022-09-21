#include "following_constraints.h"

following_constraints::following_constraints(GRBModel& md, variables var, model_type mt, parameters p, dataset& dt){
  string constraint_name;
  if (mt.base == baseModel::F){
    for (int i=0; i<p.I; i++){ // u^i_a(t),t = u^i_t,l(t) + u^i_t,r(t) + u^i_t,w
      constraint_name = "conserv_flow_N_i=" + to_string(i) + "_t=" + to_string(0);
      md.addConstr(var.u[i*p.E], GRB_EQUAL, var.u[i*p.E+1] + var.u[i*p.E+2] + var.u[i*p.E+3], constraint_name);
      for (int t=1; t<p.N; t++){
	int index_u_att = i*p.E + 3*((t-1)/2)+ 1 + ((t-1)%2),
	  index_u_tlt = i*p.E + t*3 + 1,
	  index_u_trt = i*p.E + t*3 + 2,
	  index_u_tw = i*p.E + t*3 + 3;
	constraint_name = "conserv_flow_N_i=" + to_string(i) +"_t="+to_string(t);
	md.addConstr(var.u[index_u_att], GRB_EQUAL, var.u[index_u_tlt] + var.u[index_u_trt] + var.u[index_u_tw],constraint_name);
      }
    }
    for (int i=0; i<p.I; i++){ // u^i_a(t),t = u^i_t,w
      for (int t=p.N ; t<p.N+p.L; t++){
	int index_u_tw = i*p.E + p.N*3 + t-p.N +1,
	  index_u_att = i*p.E + 3*((t-1)/2)+ 1 + ((t-1)%2);
	constraint_name = "conserv_flow_L_i=" + to_string(i) + "_t=" + to_string(t);
	md.addConstr(var.u[index_u_att], GRB_EQUAL, var.u[index_u_tw], constraint_name);
      }
    }
    for (int i=0; i<p.I; i++){ // u^i_t,w <= c_yi,t
      for (int t=0; t<p.N; t++){
	int index_u_tw = i*p.E + t*3 + 3;
	constraint_name = "uitw_LEQ_cyit_i=" + to_string(i) + "_t=" + to_string(t);
	md.addConstr(var.u[index_u_tw], GRB_LESS_EQUAL, var.c[t*p.K + dt.Y[i]], constraint_name);
      }
      for (int t = p.N; t < p.N+p.L; t++) {
	int index_u_tw = i*p.E + p.N*3 + t-p.N +1;
	constraint_name = "uitw_LEQ_cyit_i=" + to_string(i) + "_t=" + to_string(t);
	md.addConstr(var.u[index_u_tw], GRB_LESS_EQUAL, var.c[t * p.K + dt.Y[i]], constraint_name);
      }
    }
    for (int i=0; i<p.I; i++){ // u^i_t,r(t) <= sum a_j,t/d_t et u^i_t,r(t) <= sum a_j,t/d_t
      for (int t=0; t<p.N; t++){
	GRBLinExpr sum_a_jt = 0;
	if (mt.univ){
	  for (int j = 0; j < p.J; j++) {
	    sum_a_jt += var.a[t*p.J+j];
	  }
	}
	else{
	  sum_a_jt += var.d[t];
	}
	int index_u_tlt = i*p.E + t*3 + 1,
	  index_u_trt = i*p.E + t*3 + 2;
	    
	constraint_name = "uitl(t)_LEQ_sum_ajt_i=" + to_string(i) + "_t=" + to_string(t);
	md.addConstr(var.u[index_u_tlt], GRB_LESS_EQUAL, sum_a_jt, constraint_name);
	    
	constraint_name = "uitr(t)_LEQ_sum_ajt_i=" + to_string(i) + "_t=" + to_string(t);
	md.addConstr(var.u[index_u_trt], GRB_LESS_EQUAL, sum_a_jt, constraint_name);
      }
    }
  }
  else{
    for (int i=0; i<p.I; i++){ // z_it <= l_t
      for (int t=0; t<p.L; t++){
	constraint_name = "zit_LEQ_lt_i=" + to_string(i) + "_t="+to_string(t);
	md.addConstr(var.z[i*p.L+t], GRB_LESS_EQUAL,var.l[t],constraint_name);
      }
    }
    for (int t=0; t<p.L; t++){ // sum_i z_it >= Nmin*l_t
      GRBLinExpr sum_z_it = 0;
      for (int i=0; i<p.I; i++){
	if (dt.weightedPoints){
	  sum_z_it += dt.weights[i] * var.z[i*p.L+t];
	}
	else{
	  sum_z_it += var.z[i*p.L+t];
	}
      }
      constraint_name = "sumzit_GEQ_Nminlt_t=" + to_string(t);
      md.addConstr(sum_z_it, GRB_GREATER_EQUAL, p.Nmin*var.l[t], constraint_name);
    }
    for (int i=0; i<p.I; i++){ // sum_t z_it = 1
      GRBLinExpr sum_z_it = 0;
      for (int t=0; t<p.L; t++){
	sum_z_it += var.z[i*p.L+t];
      }
      constraint_name = "sumzit_EQ_1_i=" + to_string(i);
      md.addConstr(sum_z_it,GRB_EQUAL,1, constraint_name);
    }
    for (int t=0; t<p.L; t++){ // sum_k c_kt = l_t
      GRBLinExpr sum_c_kt = 0;
      for (int k=0; k<p.K; k++){
	sum_c_kt += var.c[t*p.K+k];
      }
      constraint_name = "sumckt_EQ_lt_t=" + to_string(t);
      md.addConstr(sum_c_kt, GRB_EQUAL, var.l[t], constraint_name);
    }
  }
}
