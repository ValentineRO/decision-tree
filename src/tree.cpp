#include "tree.h"

Tree::Tree(){
    D = 0;
    J = 0;
    K = 0;
}

Tree::Tree(int D_, int J_, int K_){
    D = D_;
    J = J_;
    K = K_;

    L = (int)pow(2,D);
    N = (int)pow(2,D) - 1;
    E = 1 + 3*N + L;
    this->get_tree_structure();

    a = new double[J*N];
    b = new double[N];
    eps = new double[N];

    for (int t=0; t<N; t++){
      eps[t] = 0;
    }

    c = new int[L+N];
}

Tree Tree::copy(){
  Tree newTree = Tree(D,J,K);

  for (int t=0; t<N; t++){
    for (int j=0; j<J; j++){
      newTree.a[t*J+j] = a[t*J+j];
    }
    newTree.b[t] = b[t];
    newTree.eps[t] = eps[t];
  }
  for (int t=0; t<N+L; t++){
    newTree.c[t] = c[t];
  }

  return newTree;
}

void Tree::get_tree_structure(){
  L = (int)pow(2,D);
  N = (int)pow(2,D) - 1;
  E = 1 + 3*N + L;

  /*
  A_L = new vector<int>[L];
  A_R = new vector<int>[L];

  for (int t=0; t<L; t++){
    A_L[t] = {};
    A_R[t] = {};
    int node = t+L;
    int branch;
    while (node != 1){
      branch = node % 2;
      node = node / 2;
      if(branch == 1){A_R[t].push_back(node-1);}
      else{A_L[t].push_back(node-1);}
    }
  }
  */
  leaves = new vector<int>[N];
  right_most_leaf = new int[N+L];
  
  for (int t=0; t<N; t++){
    int number_of_leaves = 2,
      node = t+1;
    node *= 2;
    while (node <= N){
      node *= 2;
      number_of_leaves *= 2;
    }
    node -= 1;
    for (int l=node; l<node+number_of_leaves; l++){
      leaves[t].push_back(l);
    }
    right_most_leaf[t] = node+number_of_leaves - 1;
  }

  for (int t=N; t<N+L; t++){
    right_most_leaf[t] = t;
  }
}

void Tree::correctSplitFunctions(){
  for (int t=0; t<N; t++){
    if (c[t] == -1){ // if there really is a split (otherwise we are going to divide by zero
      double sum_a = 0;
      for (int j=0; j<J; j++){
	sum_a += abs(a[t*J+j]);
      }
      if (abs(b[t]) > sum_a){
	cout << "abs(b) > sum_(abs(a)): this is not supposed to happen" << endl;
      }
      if (greaterThan(sum_a,0)){
	for (int j=0; j<J; j++){
	  a[t*J+j] /= sum_a;
	}
	b[t] /= sum_a;
      }
    }
  }
}

void Tree::post_processing_b(dataset& dt, bool missClassifCounts, float rho){
  int* pred = new int[dt.I];
  predict_classes(dt,pred);
  
  double* b_min = new double[N];
  double* b_max = new double[N];
  bool* is_reached = new bool[N];

  for (int t=0; t<N; t++){
    b_min[t] = 0;
    b_max[t] = 1;
    is_reached[t] = false;
  }

  for (int i=0; i<dt.I; i++){
    if (missClassifCounts||((!missClassifCounts)&&(pred[i]==dt.Y[i]))){
      int node = 0;
      for (int d=0; d<D; d++){
	if (c[node]>=0){
	  break;
	}
	is_reached[node] =true;
	
	double sum = 0;
	for (int j=0; j<J;j++){
	  sum += dt.X[i*J+j]*a[node*J+j];
	}
	if (lessThan(sum,b[node])){
	  b_min[node] = max(b_min[node],sum);
	  node = 2*node + 1;
	}
	else{
	  b_max[node] = min(b_max[node],sum);
	  node = node*2 +2;
	}
      }
    }  
  }

  for (int t=0; t<N; t++){
    if (is_reached[t]){
      b[t] = (1-rho)*b_min[t] + rho*b_max[t];
      eps[t] = min(b[t] - b_min[t], b_max[t] - b[t]);
    }
    else{
      b[t] = 0.0;
      eps[t] = 0.0;
    }
  }
    
}

void Tree::compute_ILIR(vector<int>* I_L, vector<int>* I_R, dataset& dt,  bool missClassifCounts){
  int* pred = new int[dt.I];
  predict_classes(dt,pred);
  
  for (int i=0; i<dt.I; i++){
    if (missClassifCounts||((!missClassifCounts)&&(pred[i]==dt.Y[i]))){
      int node = 0;
      for (int d=0; d<D; d++){
	if (c[node]>=0){
	  break;
	}	
	double sum = 0;
	for (int j=0; j<J;j++){
	  sum += dt.X[i*J+j]*a[node*J+j];
	}
	if (lessThan(sum,b[node])){
	  I_L[node].push_back(i); 
	  node = 2*node+1;
	}
	else{
	  I_R[node].push_back(i);
	  node = node*2 +2;
	}
      }
    }
  }
}

void Tree::post_processing_a_b(dataset& dt, bool missClassifCounts){
  int* pred = new int[dt.I];
  predict_classes(dt,pred);

  vector<int>* I_L = new vector<int>[N];
  vector<int>* I_R = new vector<int>[N];
  compute_ILIR(I_L, I_R, dt, missClassifCounts);
  
  GRBEnv env = GRBEnv();
  for (int t=0; t<N; t++){    
    if (I_L[t].size() !=0 && I_R[t].size() !=0){
      int sum_sj_star = 0;
      for (int j=0; j<J; j++){
	sum_sj_star += a[t*J+j] != 0.0;
      }
      
      GRBModel m = GRBModel(env);

      GRBVar b_var;
      GRBVar* a_var;
      GRBVar* e_var;
      GRBVar e_min;
      GRBVar* s_var;
      
      double* lb_a = new double[J];
      double* ub_a = new double[J];
      double* coef_a = new double[J];
      char* type_a = new char[J];
      string* name_a = new string[J];
      
      double* lb_e = new double[I_L[t].size()+ I_R[t].size()];
      double* ub_e = new double[I_L[t].size()+ I_R[t].size()];
      double* coef_e = new double[I_L[t].size()+ I_R[t].size()];
      char* type_e = new char[I_L[t].size()+ I_R[t].size()];
      string* name_e = new string[I_L[t].size()+ I_R[t].size()];
      
      b_var = m.addVar(-1.0, 1.0, 0.0, GRB_CONTINUOUS, "b");
      s_var = m.addVars(J, GRB_BINARY);
      for (int j=0; j<J; j++){
	lb_a[j] = -1.0;
	ub_a[j] = 1.0;
	coef_a[j] = 0.0;
	type_a[j] = GRB_CONTINUOUS;
	name_a[j] = "a" + to_string(j);
      }
      a_var = m.addVars(lb_a,ub_a,coef_a,type_a,name_a,J);
      for (int i=0; i<I_L[t].size()+ I_R[t].size(); i++){
	lb_e[i] = 0.0;
	ub_e[i] = 2.0;
	coef_e[i] = 0.0;
	type_e[i] = GRB_CONTINUOUS;
	name_e[i] = "e" + to_string(i);
      }
      e_var = m.addVars(lb_e,ub_e,coef_e,type_e,name_e,I_L[t].size()+ I_R[t].size());
      e_min = m.addVar(0.0, 2.0, 0.0, GRB_CONTINUOUS, "e_min");
      
      string constraint_name;
      constraint_name = "def_b1";
      m.addConstr(b_var,GRB_GREATER_EQUAL,-1,constraint_name);
      constraint_name = "def_b2";
      m.addConstr(b_var,GRB_LESS_EQUAL,1,constraint_name);
      
      GRBLinExpr sum_sj = 0;
      for (int j=0; j<J;j++){
	constraint_name = "a_j_geq_-sj_j="+to_string(j);
	m.addConstr(a_var[j],GRB_GREATER_EQUAL,-s_var[j],constraint_name);
	constraint_name = "aj_leq_sj_j="+to_string(j);
	m.addConstr(a_var[j],GRB_LESS_EQUAL,s_var[j],constraint_name);
	sum_sj += s_var[j];
      }
      constraint_name = "lesser_or_equal_complexity";
      m.addConstr(sum_sj,GRB_LESS_EQUAL,sum_sj_star,constraint_name);
      
      for (int i=0;i<I_L[t].size();i++){
	constraint_name = "e_min_LEQ_ei_i=" + to_string(i);
	m.addConstr(e_min,GRB_LESS_EQUAL,e_var[i],constraint_name);
	
	GRBLinExpr sum_ax = 0;
	for (int j=0; j<J; j++){
	  sum_ax += a_var[j]*dt.X[I_L[t][i]*J+j];
	}
	constraint_name = "def_ei_i="+to_string(i);
	m.addConstr(e_var[i],GRB_EQUAL,-sum_ax + b_var,constraint_name);
      }
      
      for (int i=0;i<I_R[t].size();i++){
	constraint_name = "e_min_LEQ_ei_i=" + to_string(i+I_L[t].size());
	m.addConstr(e_min,GRB_LESS_EQUAL,e_var[i+I_L[t].size()],constraint_name);
	
	GRBLinExpr sum_ax = 0;
	for (int j=0; j<J; j++){
	  sum_ax += a_var[j]*dt.X[I_R[t][i]*J+j];
	}
	constraint_name = "def_ei_i="+to_string(i+I_L[t].size());
	m.addConstr(e_var[i+I_L[t].size()],GRB_EQUAL,sum_ax - b_var,constraint_name);
      }
      GRBLinExpr expr = e_min;
      m.setObjective(expr,GRB_MAXIMIZE);

      if (setPrecision){
	m.set(GRB_DoubleParam_IntFeasTol, GLOBAL_PRECISION); // val par def : 10-5, on passe à 10-6
	m.set(GRB_DoubleParam_FeasibilityTol, GLOBAL_PRECISION); // comme la valeur par defaut
	m.set(GRB_DoubleParam_OptimalityTol, GLOBAL_PRECISION); // comme la valeur par defaut
      }
      
      //m.write("model.lp");
      //freopen("post_pr.txt", "a", stdout);
      m.optimize();
      //freopen("/dev/tty", "w", stdout);
      //remove("post_pr.txt");
      b[t] = b_var.get(GRB_DoubleAttr_X);
      eps[t] = e_min.get(GRB_DoubleAttr_X);
      double sum_a = 0;
      for (int j=0; j<J; j++){
	a[t*J+j] = a_var[j].get(GRB_DoubleAttr_X);
	sum_a += a_var[j].get(GRB_DoubleAttr_X);
      }
      if (abs(sum_a)>=1){
	b[t] /= abs(sum_a);
	for (int j=0; j<J; j++){
	  a[t*J+j] /= abs(sum_a);
	}
      }
    }
    else{
      b[t] = 0.0;
      eps[t] = 0.0;
      for (int j=0; j<J; j++){
	a[t*J+j] = 0;
      }
    }
  }
}

Tree Tree::bigger_tree(int new_D){
  Tree new_T = Tree(new_D,J,K);
    
  for (int t=0; t<N; t++){

    for (int j=0; j<J; j++){
      new_T.a[t*J+j] = a[t*J+j];
    }
    new_T.b[t] = b[t];
    new_T.eps[t] = eps[t];
  }
  for (int t=N; t<new_T.N; t++){
    for (int j=0; j<J; j++){
      new_T.a[t*J+j] = 0.0;
    }
    new_T.b[t] = 0.0;
    new_T.eps[t] = 0.0;
  }

  for (int t=0; t<N+L; t++){
    new_T.c[t] = c[t];
  }
  for (int t=N+L; t<new_T.N+new_T.L; t++){
    new_T.c[t] = -1;
  }

  return new_T;
}

Tree Tree::reduceComplexity(dataset& dt, int Nmin, double mu){
  vector<int>* pointsInNode = new vector<int>[N+L];
  this->data_points_in_last_split(dt, pointsInNode);

  Tree best_tree;
  int misclassif = dt.I;

  for (int t=0; t<N; t++){
    int left = 2*t+1,
      right = 2*t+2;
    bool is_last_split = (c[left] >= 0) && (c[right] >=0);
    if (is_last_split){
      Tree T = this->copy();
      bool changed = T.changeSplit(t, pointsInNode[t], misclassif, dt, Nmin, mu);
      if (changed){
	best_tree = T;
      }
    }
  }

  return best_tree;
}

bool Tree::changeSplit(int t, vector<int> points, int& misclassif, dataset& dt, int Nmin, double mu){
  int C = -1;
  for (int j=0; j<J; j++){
    if (a[t*J+j] !=0){
      C += 1;
    }
  }

  if (C == 0){
    int* nbPoints = new int[K];
    for (int i=0;i<points.size();i++){
      nbPoints[dt.Y[points[i]]] += 1;
    }
    int cl;
    int max = 0;
    for (int k=0; k<K; k++){
      if (nbPoints[k] > max){
	cl = k;
	max = nbPoints[k];
      }
    }
    if (misclassif > points.size() - max){
      c[2*t+1] = -1;
      c[2*t+2] = -1;
      c[t] = cl;
      for (int j=0; j<J; j++){
	a[t*J+j] = 0;
      }
      b[t] = 0;
      eps[t] = 0;
      return true;
    }
    else{
      return false;
    }
  }
  
  
  GRBEnv env = GRBEnv();
  GRBModel md = GRBModel(env);

  GRBVar* var_a = md.addVars(J, GRB_CONTINUOUS);
  GRBVar* ah = md.addVars(J, GRB_CONTINUOUS);
  GRBVar* s = md.addVars(J, GRB_BINARY);
  GRBVar var_b = md.addVar(-1.0, 1.0, 0.0, GRB_CONTINUOUS, "b");
  GRBVar* var_c = md.addVars(K*2, GRB_BINARY); // c[t*K+k]
  GRBVar d = md.addVar(0, 1, 0, GRB_BINARY, "d");
  GRBVar* l = md.addVars(2, GRB_BINARY);
  GRBVar* Lt = md.addVars(2, GRB_CONTINUOUS);
  GRBVar* Nt = md.addVars(2, GRB_CONTINUOUS);
  GRBVar* Nkt = md.addVars(2*K, GRB_CONTINUOUS); 
  GRBVar* z = md.addVars(2*points.size(), GRB_BINARY); // z[i*2+t]

  md.addConstr(-d, GRB_LESS_EQUAL, var_b, "b_GEQ_-d=");
  md.addConstr(var_b, GRB_LESS_EQUAL, d, "b_LEQ_d=");
    
  GRBLinExpr sum_ah = 0,
    sum_s = 0;
  
  for (int j=0; j<J; j++){
    md.addConstr(-ah[j], GRB_LESS_EQUAL, var_a[j], "-ah_LEQ_a_j="+to_string(j));
    md.addConstr(var_a[j], GRB_LESS_EQUAL, ah[j], "a_LEQ_ah_j="+to_string(j));
    md.addConstr(-s[j], GRB_LESS_EQUAL, var_a[j], "-s_LEQ_a_j="+to_string(j));
    md.addConstr(var_a[j], GRB_LESS_EQUAL, s[j], "a_LEQ_s_j="+to_string(j));
    md.addConstr(s[j], GRB_LESS_EQUAL, d, "s_LEQ_d_j="+to_string(j));
    sum_ah += ah[j];
    sum_s += s[j];
  }
  md.addConstr(sum_ah,GRB_LESS_EQUAL,d, "sum_ah_LEQ_d");
  md.addConstr(sum_s,GRB_GREATER_EQUAL,d, "sum_s_GEQ_d");
  md.addConstr(sum_s,GRB_LESS_EQUAL,C, "sum_s_LEQ_C");

  for (int lf=0; lf<2; lf++){
    GRBLinExpr sum_z = 0;
    for (int i=0; i<points.size(); i++){
      md.addConstr(z[i*2+lf], GRB_LESS_EQUAL, l[lf], "z_LEQ_l_i="+to_string(points[i])+"_t="+to_string(lf));
      sum_z += z[i*2+lf];
    }
    md.addConstr(sum_z,GRB_GREATER_EQUAL,Nmin*l[lf], "sum_z_GEQ_Nminxl_t="+to_string(lf));
    md.addConstr(sum_z,GRB_EQUAL, Nt[lf], "def_Nt_t="+to_string(lf));
    GRBLinExpr sum_c = 0;
    for (int k=0; k<K; k++){
      sum_c += var_c[lf*K+k];
      GRBLinExpr sum_zk = 0;
      for (int i=0; i<points.size(); i++){
	if (dt.Y[points[i]] ==k){
	  sum_zk += z[i*2+lf];
	}
      }
      md.addConstr(sum_zk, GRB_EQUAL, Nkt[lf*K+k], "def_Nkt_t="+to_string(lf)+"_k="+to_string(k));
      md.addConstr(Lt[lf], GRB_GREATER_EQUAL, Nt[lf] - Nkt[lf*K+k] - dt.I*(1-var_c[lf*K+k]), "def_Lt1_t="+to_string(lf)+"_k="+to_string(k));
      md.addConstr(Lt[lf], GRB_LESS_EQUAL, Nt[lf] - Nkt[lf*K+k] + dt.I*var_c[lf*K+k], "def_Lt2_t="+to_string(lf)+"_k="+to_string(k));
    }
    md.addConstr(sum_c,GRB_EQUAL,l[lf],"sum_c_EQ_l_t="+to_string(lf));
    md.addConstr(Lt[lf], GRB_GREATER_EQUAL, 0, "def_Lt3_t="+to_string(lf));
  }

  for (int i=0; i<points.size(); i++){
    md.addConstr(z[i*2]+z[i*2+1], GRB_EQUAL, 1, "sum_z_EQ_1_i="+to_string(points[i]));
    GRBLinExpr sum_ax = 0;
    for (int j=0; j<J; j++){
      sum_ax += var_a[j]*dt.X[points[i]*J+j];
    }
    md.addConstr(sum_ax+mu,GRB_LESS_EQUAL, var_b + (2+mu)*(1-z[i*2]), "branching_left_i="+to_string(points[i]));
    md.addConstr(sum_ax,GRB_GREATER_EQUAL, var_b - 2*(1-z[i*2+1]), "branching_right_i="+to_string(points[i]));
  }

  md.setObjective(Lt[0] + Lt[1], GRB_MINIMIZE);
  if (setPrecision){
    md.set(GRB_DoubleParam_IntFeasTol, GLOBAL_PRECISION); // val par def : 10-5, on passe à 10-6
    md.set(GRB_DoubleParam_FeasibilityTol, GLOBAL_PRECISION); // comme la valeur par defaut
    md.set(GRB_DoubleParam_OptimalityTol, GLOBAL_PRECISION); // comme la valeur par defaut
  }
  md.set("TimeLimit", to_string(60)); // just in case it takes too much time
  //freopen("reducingComplexity.txt", "a", stdout);
  md.optimize();
  //freopen("/dev/tty", "w", stdout);
  //remove("reducingComplexity.txt");

  if (misclassif > md.get(GRB_DoubleAttr_ObjVal)){
    misclassif = md.get(GRB_DoubleAttr_ObjVal);
    b[t] = var_b.get(GRB_DoubleAttr_X);
    eps[t] = 0.0;
    for (int j=0; j<J; j++){
      a[t*J+j] = var_a[j].get(GRB_DoubleAttr_X);
    }
    for (int k=0; k<K; k++){
      if (var_c[k].get(GRB_DoubleAttr_X) >= 0.99){
	c[2*t+1] = k;
      }
      if (var_c[K+k].get(GRB_DoubleAttr_X) >= 0.99){
	c[2*t+2] = k;
      }
    }
    return true;
  }
  return false;
}

Tree Tree::removeSplit(dataset& dt){
  int best_tree_ind = -1;
  int class_of_new_leaf = -1;
  int misclassif = dt.I;

  int* rep = new int[(N+L)*K];
  this->data_points_per_leaves(dt, rep);
  
  for (int t=0; t<N; t++){
    // let's identify is the node t is a final split ie if t has 2 leaves
    int left = 2*t+1,
      right = 2*t+2;
    bool is_last_split = (c[left] >= 0) && (c[right] >=0);
    if (is_last_split){
      int nbPoints = 0,
	nbMaxPoints = 0,
	cl = -1;
      for (int k=0; k<K; k++){
	nbPoints += rep[t*K+k];
	if (nbMaxPoints < rep[t*K+k]){
	  cl = k;
	  nbMaxPoints = rep[t*K+k];
	}
      }
      if (nbPoints - nbMaxPoints < misclassif){
	misclassif = nbPoints - nbMaxPoints;
	class_of_new_leaf = cl;
	best_tree_ind = t;
      }
    }
  }
  Tree best_tree = this->copy();
  
  best_tree.c[2*best_tree_ind+1] = -1;
  best_tree.c[2*best_tree_ind+2] = -1;
  best_tree.c[best_tree_ind] = class_of_new_leaf;
  for (int j=0; j<J; j++){
    best_tree.a[best_tree_ind*J+j] = 0;
  }
  best_tree.b[best_tree_ind] = 0;
  best_tree.eps[best_tree_ind] = 0;
  
  return best_tree;
}

int Tree::predict_class(double x[]){
  int node = 1;
  for (int d=0; d<D;d++){
    if (c[node-1] >= 0){
      return c[node-1];
    }
    double sum=0;
    for (int j=0; j<J;j++){
      sum += x[j]*a[(node-1)*J+j];
    }
    if (lessThan(sum,b[node-1])){node = 2*node;}
    else{node = node*2 +1;}
    
  }
  return c[node-1];
}

void Tree::predict_classes(dataset& dt,int predictions[]){
  for (int i=0;i<dt.I;i++){
    int node = 0;
    for (int d=0; d<D;d++){
      if (c[node]>=0){
	break;
      }
      double sum=0;
      for (int j=0; j<J;j++){
	sum += dt.X[i*J+j]*a[node*J+j];
      }
      if (lessThan(sum,b[node])){
	node = 2*node+1;
      }
      else{
	node = 2*node+2;
      }
    }
    predictions[i] = c[node];
  }
}
int Tree::prediction_errors(dataset& dt){
  int* pred = new int[dt.I];
  predict_classes(dt,pred);

  int errors = 0;
  for (int i=0; i<dt.I; i++){
    if (dt.Y[i]!=pred[i]){
      if (dt.weightedPoints){
	errors += dt.weights[i];
      }
      else{
	errors += 1;	
      }
    }
  }

  return errors;
}

void Tree::data_points_per_leaves(dataset& dt, int repartition[]){
  for (int t=0; t<N+L; t++){
    for (int k=0; k<K; k++){
      repartition[t*K+k] = 0;
    }
  }
  
  for (int i=0; i<dt.I; i++){
    int node = 0;
    for (int d=0; d<D; d++){
      if (c[node-1]>=0){
	break;
      }
      double sum=0;
      
      for (int j=0; j<J;j++){
	sum += dt.X[i*J+j]*a[node*J+j];
	
      }
      
      if (lessThan(sum,b[node])){
	node = 2*node+1;
      }
      else{
	node = 2*node+2;
      }
      
    }
    repartition[node*K+dt.Y[i]] += 1;
  }
}

void Tree::data_points_in_last_split(dataset& dt, vector<int> points_in_leaf[]){
  for (int i=0; i<dt.I; i++){
    int node = 0;
    for (int d=0; d<D; d++){
      if (c[node]>=0){
	break;
      }
      double sum=0;
      
      for (int j=0; j<J;j++){
	sum += dt.X[i*J+j]*a[node*J+j];
	
      }
      
      if (lessThan(sum,b[node])){
	node = 2*node+1;
      }
      else{
	node = node*2 +2;
      }
      
    }
    points_in_leaf[(node-1)/2].push_back(i);
  }
  
}

void Tree::predict_leaves(dataset& dt, int leaves[]){
  for (int i=0;i<dt.I;i++){
    int node = 1;
    for (int d=0; d<D;d++){
      if (c[node-1]>=0){
	break;
      }
      double sum=0;
      for (int j=0; j<J;j++){
	sum += dt.X[i*J+j]*a[(node-1)*J+j];
	
      }
      if (lessThan(sum,b[node-1])){node = 2*node;}
      else{node = node*2 +1;}
    }
    if (node-1<N){
      leaves[i] = right_most_leaf[node-1]-N;
    }
    else{
      leaves[i] = node-1-N;
    }
  }
}

void Tree::predict_paths(dataset& dt, vector<int> paths[]){
  for (int i=0; i<dt.I; i++){
    int node = 1;
    paths[i].push_back(0);
    for (int d=0; d<D;d++){
      if (c[node-1]>=0){
	break;
      }
      double sum=0;
      for (int j=0; j<J;j++){
	sum += dt.X[i*J+j]*a[(node-1)*J+j];
      }
      if (lessThan(sum,b[node-1])){node = 2*node;}
      else{node = 2*node+1;}
      paths[i].push_back(node-1);
    }
  }
}

void Tree::write_tree(string namefile){
    fstream file;
    file.open(namefile,ios::out);
    file << D << "\n";
    file << J << "\n";
    file << K << "\n";

    for (int t =0;t<N;t++){
        for (int j=0;j<J;j++){
            file << a[t*J+j]<<"\n";
        }
        file <<b[t]<<"\n";
	file << eps[t] << "\n";
	file << c[t] << "\n";
        
    }
    for (int t=0;t<L;t++){
        file << c[t+N]<< "\n";
    }
    file.close();
}

Tree::Tree(string namefile){
  fstream file;
  file.open(namefile,ios::in);

  file >> D;
  file >> J;
  file >> K;
    
  L = (int)pow(2,D);
  N = (int)pow(2,D) - 1;
  E = 1 + 3*N + L;
  this->get_tree_structure();

  a = new double[J*N];
  b = new double[N];
  eps = new double[N];

  c = new int[L+N];

  for (int t=0; t<N; t++){
    for (int j=0; j<J; j++){
      file >> a[t*J+j];
    }
    file >> b[t];
    file >> eps[t];
    file >> c[t];
  }

  for (int t=0; t<L; t++){
    file >> c[t+N];
  }
  file.close();
}
