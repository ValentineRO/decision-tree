#include "counting_errors.h"

counting_errors::counting_errors(GRBModel& md, variables var, model_type mt, parameters p, dataset& dt){
  string constraint_name;
  if (mt.base==baseModel::OCT){
    for (int t=0; t<p.L; t++){ // N_kt = sum Y_ik*z_it
      for (int k=0; k<p.K; k++){
	GRBLinExpr sum_zit = 0;
	for (int i=0; i<p.I; i++){
	  if (dt.Y[i] == k){
	    if (dt.weightedPoints){
	      sum_zit += dt.weights[i] * var.z[i*p.L+t];
	    }
	    else{
	      sum_zit += var.z[i*p.L+t];
	    }
	  }
	}
	constraint_name = "Nkt_def_t=" + to_string(t) + "_k=" + to_string(k);
	md.addConstr(var.Nkt[t*p.K+k], GRB_EQUAL, sum_zit,constraint_name);
      }
    }
    for (int t=0; t<p.L; t++){ // Nt = sum z_it
      GRBLinExpr sum_zit = 0;
      for (int i=0; i<p.I; i++){
	if (dt.weightedPoints){
	  sum_zit += dt.weights[i] * var.z[i*p.L+t];
	}
	else{
	  sum_zit += var.z[i*p.L+t];
	}
      }
      constraint_name = "Nt_def_t=" + to_string(t) ;
      md.addConstr(var.Nt[t], GRB_EQUAL, sum_zit, constraint_name);
    }
    for (int t=0; t<p.L; t++){ // L_t >= Nt - N_kt - I*(1-c_kt)
      for (int k=0; k<p.K; k++){
	constraint_name = "Lt_def1_t=" + to_string(t) + "_k=" + to_string(k);
	if (dt.weightedPoints){
	  md.addConstr(var.Lt[t], GRB_GREATER_EQUAL, var.Nt[t] - var.Nkt[t*p.K+k] - dt.initialI*(1-var.c[t*p.K+k]), constraint_name);
	}
	else{
	  md.addConstr(var.Lt[t], GRB_GREATER_EQUAL, var.Nt[t] - var.Nkt[t*p.K+k] - p.I*(1-var.c[t*p.K+k]), constraint_name);
	}
      }
    }
    for (int t=0; t<p.L; t++){ // L_t <= Nt - N_kt + I*c_kt
      for (int k=0; k<p.K; k++){
	constraint_name = "Lt_def2_t=" + to_string(t) + "_k=" + to_string(k);
	if (dt.weightedPoints){
	  md.addConstr(var.Lt[t], GRB_LESS_EQUAL, var.Nt[t] - var.Nkt[t*p.K+k] + dt.initialI*var.c[t*p.K+k], constraint_name);
	}
	else{
	  md.addConstr(var.Lt[t], GRB_LESS_EQUAL, var.Nt[t] - var.Nkt[t*p.K+k] + p.I*var.c[t*p.K+k], constraint_name);
	}
      }
    }
    for (int t=0; t<p.L; t++){ // L_t >= 0
      constraint_name = "Lt_def3_t=" + to_string(t);
      md.addConstr(var.Lt[t], GRB_GREATER_EQUAL, 0, constraint_name);
    }
  }
  if (mt.base==baseModel::FOCT){
    for (int t=0; t<p.L; t++){
      for (int i=0; i<p.I; i++){
	for (int k=0; k<p.K; k++){
	  // theta_tik >= 0
	  constraint_name = "theta_tik_def1_t="+to_string(t)+"_i="+to_string(i)+"_k="+to_string(k);
	  md.addConstr(var.theta[t*p.I*p.K+i*p.K+k],GRB_GREATER_EQUAL,0,constraint_name);

	  // theta_tik >= c_tk + z_it - 1
	  constraint_name = "theta_tik_def2_t="+to_string(t)+"_i="+to_string(i)+"_k="+to_string(k);
	  md.addConstr(var.theta[t*p.I*p.K+i*p.K+k],GRB_GREATER_EQUAL,var.c[t*p.K+k]+var.z[i*p.L+t]-1,constraint_name);
	}
      }
    }
  }
  if (mt.base==baseModel::GOCT){
    for (int t=0; t<p.L; t++){
      for (int k=0; k<p.K; k++){
	GRBLinExpr sum_zit = 0;
	int I_Ik = 0;
	for (int i=0; i<p.I; i++){
	  if (dt.Y[i] != k){
	    if (dt.weightedPoints){
	      I_Ik += dt.weights[i];
	      sum_zit += dt.weights[i] * var.z[i*p.L+t];
	    }
	    else{
	      I_Ik += 1;
	      sum_zit += var.z[i*p.L+t];
	    }
	  }
	}
	// theta_tk >= 0
	constraint_name = "theta_tk_def1_t="+to_string(t)+"_k="+to_string(k);
	md.addConstr(var.theta[t*p.K+k],GRB_GREATER_EQUAL,0,constraint_name);

	// theta_tk >= sum_zit - |I \ Ik|(1- c_kt)
	constraint_name = "theta_tk_def2_t="+to_string(t)+"_k="+to_string(k);
	md.addConstr(var.theta[t*p.K+k],GRB_GREATER_EQUAL, sum_zit - I_Ik * var.c[t*p.K+k],constraint_name);
      }
    }
  }
  
}
