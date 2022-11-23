#include "objective.h"

objective::objective(GRBModel& md, variables var, model_type mt, parameters p, dataset& dt){
  // computation of the complexity of the tree:
  GRBLinExpr complx = 0;
  if (mt.univ){
    if (mt.base == baseModel::F){
      for (int j=0; j<p.J; j++){
	for (int t=0; t<p.N; t++){
	  complx += var.a[t*p.J+j];
	}
      }
    }
    else{
      for (int t=0; t<p.N; t++){
	complx += var.d[t];
      }
    }
  }
  else{
    for (int j=0; j<p.J; j++){
      for (int t=0; t<p.N; t++){
	complx += var.s[t*p.J+j];
      }
    }
  }

  if (mt.C_const){ // penalizing the complexity with a constraint
    md.addConstr(complx,GRB_LESS_EQUAL,p.C,"limit_complexity");
  }

  // we have to differentiate between the quad obj and the rest (lin obj)
  if (mt.base == baseModel::QOCT){
    GRBQuadExpr expr = 0;
    for (int t = 0; t < p.L; t++) {
      for (int k = 0; k < p.K; k++) {
	for (int i = 0; i < dt.I; i++) {
	  if (dt.Y[i] != k) {
	    if (dt.weightedPoints){
	      expr += dt.weights[i] * var.c[t * p.K + k] * var.z[i * p.L + t] / p.L_hat;
	    }
	    else{
	      expr += var.c[t * p.K + k] * var.z[i * p.L + t] / p.L_hat;
	    }
	  }
	}
      }
    }

    expr += complx*p.alph;

    if (mt.eps){
      double xi = 1/((1+(int)(!mt.univ))*p.L_hat*p.N);
      for (int t=0; t<p.N; t++){
	expr += xi*var.eps[t];
      }
    }
    md.setObjective(expr,GRB_MINIMIZE);
  }
  else{
    GRBLinExpr expr = 0;
    if (mt.base == baseModel::OCT){
      for (int t=0; t<p.L; t++){
        expr += var.Lt[t]/ p.L_hat;
      }
    }
    if (mt.base == baseModel::F){
      if (dt.weightedPoints){
	expr += (float)dt.initialI/ (float)p.L_hat;
      }
      else{
	expr += (float)p.I/ (float)p.L_hat;
      }
      for (int i=0; i<p.I; i++){
	if (dt.weightedPoints){
	  expr -= dt.weights[i] * var.u[i*p.E]/ p.L_hat;
	}
	else{
	  expr -= var.u[i*p.E]/ p.L_hat;
	}
      }
    }
    if (mt.base == baseModel::FOCT){
      for (int i=0; i<p.I; i++){
	for (int k=0; k<p.K; k++){
	  if (dt.Y[i] !=k){
	    for (int t=0; t<p.L; t++){
	      if (dt.weightedPoints){
		expr += dt.weights[i] * var.theta[t*p.I*p.K+i*p.K+k]/ p.L_hat;
	      }
	      else{
		expr += var.theta[t*p.I*p.K+i*p.K+k]/ p.L_hat;
	      }
	    }
	  }
	}
      }
    }
    if (mt.base == baseModel::GOCT){
      for (int k=0; k<p.K; k++){
	for (int t=0; t<p.L; t++){
	  expr += var.theta[t*p.K+k]/ p.L_hat;
	}
      }
    }
    
    expr += complx*p.alph;
    
    if (mt.eps){
      double xi = 1/((1+(int)(!mt.univ))*p.L_hat*p.N);
      for (int t=0; t<p.N; t++){
	expr += xi*var.eps[t];
      }
    }
    md.setObjective(expr,GRB_MINIMIZE);
  }
}
