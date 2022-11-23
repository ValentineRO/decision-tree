#include "models.h"

build_model::build_model(GRBModel& md, dataset& dt, model_type modelt, parameters p){  
  mt = modelt;
  param = p;
  var = variables(md, mt, param);

  md.update();

  obj = objective(md, var, mt, param, dt);
  flwc = following_constraints(md, var, mt, param, dt);
  tsc = tree_structuring_constraints(md, var, mt, param);
  bc = branching_constraints(md, var, mt, param, dt);
  ce = counting_errors(md , var, mt, param, dt);

  if (mt.base==baseModel::QOCT){
    md.set(GRB_IntParam_NonConvex,2);
    //    md.set(GRB_IntParam_PreQLinearize,1);
  }

  md.update();

  if (setPrecision){
    md.set(GRB_DoubleParam_IntFeasTol, GLOBAL_PRECISION); // val par def : 10-5, on passe à 10-6
    md.set(GRB_DoubleParam_FeasibilityTol, GLOBAL_PRECISION); // comme la valeur par defaut
    md.set(GRB_DoubleParam_OptimalityTol, GLOBAL_PRECISION); // comme la valeur par defaut
  }
  
  //md.write("model.lp");
}

void build_model::add_warmstart(Tree T, dataset &dt){
  if (T.D < param.D){
    T = T.bigger_tree(param.D);
  }
  int* val_d = new int[param.N];
  for (int t=0; t<param.N; t++){
    var.b[t].set(GRB_DoubleAttr_Start,T.b[t]);
    double sum_a = 0;
    for (int j=0; j<T.J; j++){
      var.a[t*T.J+j].set(GRB_DoubleAttr_Start,T.a[t*T.J+j]);
      sum_a += abs(T.a[t*T.J+j]);
    }
    val_d[t] = (int)(sum_a>0); // dt is also useful for Funiv but not defined in the model
    if (mt.base != baseModel::F or !mt.univ){
      if (sum_a>0){
	var.d[t].set(GRB_DoubleAttr_Start,1);
      }
      else{
	var.d[t].set(GRB_DoubleAttr_Start,0);
      }
    }
    if (!mt.univ){
      for (int j=0; j<T.J; j++){
	if (T.a[t*T.J+j] !=0){
	  var.a_h[t*T.J+j].set(GRB_DoubleAttr_Start,abs(T.a[t*T.J+j]));
	  var.s[t*T.J+j].set(GRB_DoubleAttr_Start,1);
	}
	else{
	  var.a_h[t*T.J+j].set(GRB_DoubleAttr_Start,0);
	  var.s[t*T.J+j].set(GRB_DoubleAttr_Start,0);
	}
      }
    }
  }
  if (mt.base == baseModel::F){
    for (int t=0; t<param.L+param.N; t++){
      if (T.c[t]>=0){
	for (int k=0; k<T.K; k++){
	  var.c[t*T.K+k].set(GRB_DoubleAttr_Start,(int)(T.c[t]==k));
	}
      }
      else{
	if (t>=param.N){ //default value is class 0
	  var.c[t*T.K].set(GRB_DoubleAttr_Start,1);
	}
	else{
	  if (val_d[t] >= 0.99){
	    var.c[t*T.K].set(GRB_DoubleAttr_Start,0);
	  }
	  else{
	    var.c[t*T.K].set(GRB_DoubleAttr_Start,1);
	  }
	}
      	for (int k=1; k<T.K; k++){
	  var.c[t*T.K+k].set(GRB_DoubleAttr_Start,0);
	}
      }
    }
    vector<int>* paths = new vector<int>[dt.I];
    T.predict_paths(dt,paths);

    int* preds = new int[dt.I];
    if (mt.univ){
      T.predict_classes(dt,preds);
    }
    else{
      T.predict_classes(dt,preds);
    }
    
    for (int i=0; i<dt.I; i++){
      if (preds[i] == dt.Y[i]){
	int* arcs_in_path = new int[param.E];
	arcs_in_path[0] = 1;
	for (int e=1; e<param.E; e++){ // initialisation
	  arcs_in_path[e] = 0;
	}
	for (int e=0; e< paths[i].size() - 1; e++){
	  if (paths[i][e+1] % 2 ==0){ // fils droit
	    arcs_in_path[3*paths[i][e]+2] = 1;
	  }
	  else{// fils gauche
	    arcs_in_path[3*paths[i][e]+1] = 1;
	  }
	}
	if (paths[i][paths[i].size()-1] >= param.N){
	  arcs_in_path[2*param.N + paths[i][paths[i].size()-1]+1] = 1;
	}
	else{
	  arcs_in_path[3*paths[i][paths[i].size()-1]+3] = 1;
	}
	
	for (int e=0; e<param.E; e++){
	  var.u[i*param.E+e].set(GRB_DoubleAttr_Start,arcs_in_path[e]);
	}
      }
    }
  }
  else{
    int* val_c = new int[param.L];
    for (int t=0; t<param.L; t++){
      val_c[t] = T.c[t+param.N];
    }
    for (int t=0; t<param.N; t++){
      if (T.c[t] != -1){
	val_c[T.right_most_leaf[t]-param.N] = T.c[t];
      }
    }
    for (int t=0; t<param.L; t++){
      for (int k=0; k<T.K; k++){
	var.c[t*T.K+k].set(GRB_DoubleAttr_Start,(int)(val_c[t]==k)); // if there is a true leaf, then we will have a ckt=1, otherwise they all equal 0
      }
    }
    
    int* leaves = new int[param.I];
    if (mt.univ){
      T.predict_leaves(dt,leaves);
    }
    else{
      T.predict_leaves(dt,leaves);
    }
    
    for (int t=0; t<param.L; t++){
      bool is_open = false;
      for (int i=0; i<dt.I; i++){
	if (t == leaves[i]){
	  var.z[i*param.L+t].set(GRB_DoubleAttr_Start,1);
	  is_open = true;
	}
	else{
	  var.z[i*param.L+t].set(GRB_DoubleAttr_Start,0);
	}
      }
      if (is_open){
	var.l[t].set(GRB_DoubleAttr_Start,1);
      }
      else{
	var.l[t].set(GRB_DoubleAttr_Start,0);
      }
    }

    if (mt.base == baseModel::FOCT){
      for (int t=0; t<param.L; t++){
	for (int i=0; i<dt.I; i++){
	  for (int k=0; k<dt.K; k++){
	    
	    if (leaves[i]==t and val_c[t]==k){
	      var.theta[t*dt.I*dt.K+i*dt.K+k].set(GRB_DoubleAttr_Start,1);
	    }
	    else{
	      var.theta[t*dt.I*dt.K+i*dt.K+k].set(GRB_DoubleAttr_Start,0);
	    }
	  }
	}
      }
    }
    if (mt.base == baseModel::GOCT){
      for (int t=0; t<param.L; t++){
	for (int k=0; k<dt.K; k++){
	  int sum = 0;
	  for (int i=0; i<dt.I; i++){
	    if (leaves[i]==t and val_c[t]==k and dt.Y[i] != k){
	      sum += 1;
	    }
	  }
	  var.theta[t*dt.K+k].set(GRB_DoubleAttr_Start,sum);
	}
      }
    }
    if (mt.base == baseModel::OCT){
      for (int t=0; t<param.L; t++){
	int nt = 0;
	int lt = 0;
	int* nkt = new int[dt.K];
	for (int k=0; k<dt.K; k++){
	  nkt[k] = 0;
	}
	for (int i=0; i<dt.I; i++){
	  nt += (int)(leaves[i] == t);
	  lt += (int)(leaves[i] == t and val_c[t] != dt.Y[i]);
	  for (int k=0; k<dt.K; k++){
	    nkt[k] += (int)(leaves[i] == t and k == dt.Y[i]);
	  }
	}
	var.Lt[t].set(GRB_DoubleAttr_Start,lt);
	var.Nt[t].set(GRB_DoubleAttr_Start,nt);
	for (int k=0; k<dt.K; k++){
	  var.Nkt[t*dt.K+k].set(GRB_DoubleAttr_Start,nkt[k]);
	}
      }
    }
  }
}

void build_model::addWarmStartInConstraints(GRBModel& md, dataset& dt, Tree T){
  if (T.D < param.D){
    T = T.bigger_tree(param.D);
  }
  string constraint_name;
  int* val_d = new int[param.N];
  for (int t=0; t<param.N; t++){
    constraint_name = "set_var_bt_t="+to_string(t);
    md.addConstr(var.b[t],GRB_EQUAL, T.b[t], constraint_name);
    double sum_a = 0;
    for (int j=0; j<T.J; j++){
      constraint_name = "set_var_ajt_j="+to_string(j)+"_t="+to_string(t);
      md.addConstr(var.a[t*T.J+j],GRB_EQUAL, T.a[t*T.J+j], constraint_name);
      sum_a += abs(T.a[t*T.J+j]);
    }
    val_d[t] = (int)(sum_a>0); // dt is also useful for Funiv but not defined in the model
    if (mt.base != baseModel::F or !mt.univ){
      if (sum_a>0){
	constraint_name = "set_var_dt_t="+to_string(t);
	md.addConstr(var.d[t],GRB_EQUAL, 1, constraint_name);
      }
      else{
	constraint_name = "set_var_dt_t="+to_string(t);
	md.addConstr(var.d[t],GRB_EQUAL, 0, constraint_name);
      }
    }
    if (!mt.univ){
      for (int j=0; j<T.J; j++){
	if (T.a[t*T.J+j] !=0){
	  constraint_name = "set_var_ahjt_j="+to_string(j)+"_t="+to_string(t);
	  md.addConstr(var.a_h[t*T.J+j],GRB_EQUAL, abs(T.a[t*T.J+j]), constraint_name);
	  constraint_name = "set_var_sjt_j="+to_string(j)+"_t="+to_string(t);
	  md.addConstr(var.s[t*T.J+j],GRB_EQUAL, 1, constraint_name);
	}
	else{
	  constraint_name = "set_var_ahjt_j="+to_string(j)+"_t="+to_string(t);
	  md.addConstr(var.a_h[t*T.J+j],GRB_EQUAL, 0, constraint_name);
	  constraint_name = "set_var_sjt_j="+to_string(j)+"_t="+to_string(t);
	  md.addConstr(var.s[t*T.J+j],GRB_EQUAL, 0, constraint_name);
	}
      }
    }
  }
  if (mt.base == baseModel::F){
    for (int t=0; t<param.L+param.N; t++){
      if (T.c[t]>=0){
	for (int k=0; k<T.K; k++){
	  constraint_name = "set_var_ckt_k="+to_string(k)+"_t="+to_string(t);
	  md.addConstr(var.c[t*T.K+k],GRB_EQUAL, (int)(T.c[t]==k), constraint_name);
	}
      }
      else{
	if (t>=param.N){ //default value is class 0
	  constraint_name = "set_var_ckt_k="+to_string(0)+"_t="+to_string(t);
	  md.addConstr(var.c[t*T.K],GRB_EQUAL, 1, constraint_name);
	}
	else{
	  if (val_d[t] >= 0.99){
	    constraint_name = "set_var_ckt_k="+to_string(0)+"_t="+to_string(t);
	    md.addConstr(var.c[t*T.K],GRB_EQUAL, 0, constraint_name);
	  }
	  else{
	    constraint_name = "set_var_ckt_k="+to_string(0)+"_t="+to_string(t);
	    md.addConstr(var.c[t*T.K],GRB_EQUAL, 1, constraint_name);
	  }
	}
      	for (int k=1; k<T.K; k++){
	  constraint_name = "set_var_ckt_k="+to_string(k)+"_t="+to_string(t);
	  md.addConstr(var.c[t*T.K+k],GRB_EQUAL, 0, constraint_name);
	}
      }
    }
    vector<int>* paths = new vector<int>[dt.I];
    T.predict_paths(dt,paths);

    int* preds = new int[dt.I];
    if (mt.univ){
      T.predict_classes(dt,preds);
    }
    else{
      T.predict_classes(dt,preds);
    }
    
    for (int i=0; i<dt.I; i++){
      if (preds[i] == dt.Y[i]){
	int* arcs_in_path = new int[param.E];
	arcs_in_path[0] = 1;
	for (int e=1; e<param.E; e++){ // initialisation
	  arcs_in_path[e] = 0;
	}
	for (int e=0; e< paths[i].size() - 1; e++){
	  if (paths[i][e+1] % 2 ==0){ // fils droit
	    arcs_in_path[3*paths[i][e]+2] = 1;
	  }
	  else{// fils gauche
	    arcs_in_path[3*paths[i][e]+1] = 1;
	  }
	}
	if (paths[i][paths[i].size()-1] >= param.N){
	  arcs_in_path[2*param.N + paths[i][paths[i].size()-1]+1] = 1;
	}
	else{
	  arcs_in_path[3*paths[i][paths[i].size()-1]+3] = 1;
	}
	
	for (int e=0; e<param.E; e++){
	  constraint_name = "set_var_u_ie_i="+to_string(i)+"_e="+to_string(e);
	  md.addConstr(var.u[i*param.E+e],GRB_EQUAL, arcs_in_path[e], constraint_name);
	}
      }
    }
  }
  else{
    int* val_c = new int[param.L];
    for (int t=0; t<param.L; t++){
      val_c[t] = T.c[t+param.N];
    }
    for (int t=0; t<param.N; t++){
      if (T.c[t] != -1){
	val_c[T.right_most_leaf[t]-param.N] = T.c[t];
      }
    }
    for (int t=0; t<param.L; t++){
      for (int k=0; k<T.K; k++){
	constraint_name = "set_var_ckt_k="+to_string(k)+"_t="+to_string(t);
	md.addConstr(var.c[t*T.K+k],GRB_EQUAL, (int)(val_c[t]==k), constraint_name);// if there is a true leaf, then we will have a ckt=1, otherwise they all equal 0
      }
    }
    
    int* leaves = new int[param.I];
    if (mt.univ){
      T.predict_leaves(dt,leaves);
    }
    else{
      T.predict_leaves(dt,leaves);
    }
    
    for (int t=0; t<param.L; t++){
      bool is_open = false;
      for (int i=0; i<dt.I; i++){
	if (t == leaves[i]){
	  constraint_name = "set_var_z_it_i="+to_string(i)+"_t="+to_string(t);
	  md.addConstr(var.z[i*param.L+t],GRB_EQUAL, 1, constraint_name);
	  is_open = true;
	}
	else{
	  constraint_name = "set_var_z_it_i="+to_string(i)+"_t="+to_string(t);
	  md.addConstr(var.z[i*param.L+t],GRB_EQUAL, 0, constraint_name);
	}
      }
      if (is_open){
	constraint_name = "set_var_lt_t="+to_string(t);
	md.addConstr(var.l[t],GRB_EQUAL, 1, constraint_name);
      }
      else{
	constraint_name = "set_var_lt_t="+to_string(t);
	md.addConstr(var.l[t],GRB_EQUAL, 0, constraint_name);
      }
    }

    if (mt.base == baseModel::FOCT){
      for (int t=0; t<param.L; t++){
	for (int i=0; i<dt.I; i++){
	  for (int k=0; k<dt.K; k++){
	    if (leaves[i]==t and val_c[t]==k){
	      constraint_name = "set_theta_itk_i="+to_string(i)+"_t="+to_string(t)+"_k="+to_string(k);
	      md.addConstr(var.theta[t*dt.I*dt.K+i*dt.K+k],GRB_EQUAL,1,constraint_name);
	    }
	    else{
	      constraint_name = "set_theta_itk_i="+to_string(i)+"_t="+to_string(t)+"_k="+to_string(k);
	      md.addConstr(var.theta[t*dt.I*dt.K+i*dt.K+k],GRB_EQUAL,0,constraint_name);
	    }
	  }
	}
      }
    }
    if (mt.base == baseModel::GOCT){
      for (int t=0; t<param.L; t++){
	for (int k=0; k<dt.K; k++){
	  int sum = 0;
	  for (int i=0; i<dt.I; i++){
	    if (leaves[i]==t and val_c[t]==k and dt.Y[i] != k){
	      sum += 1;
	    }
	  }
	  constraint_name = "set_theta_tk_t="+to_string(t)+"_k="+to_string(k);
	  md.addConstr(var.theta[t*dt.K+k],GRB_EQUAL,sum,constraint_name);
	}
      }
    }
    if (mt.base == baseModel::OCT){
      for (int t=0; t<param.L; t++){
	int nt = 0;
	int lt = 0;
	int* nkt = new int[dt.K];
	for (int k=0; k<dt.K; k++){
	  nkt[k] = 0;
	}
	for (int i=0; i<dt.I; i++){
	  nt += (int)(leaves[i] == t);
	  lt += (int)(leaves[i] == t and val_c[t] != dt.Y[i]);
	  for (int k=0; k<dt.K; k++){
	    nkt[k] += (int)(leaves[i] == t and k == dt.Y[i]);
	  }
	}
	constraint_name = "set_var_Lt_t="+to_string(t);
	md.addConstr(var.Lt[t],GRB_EQUAL, lt, constraint_name);
	constraint_name = "set_var_Nt_t="+to_string(t);
	md.addConstr(var.Nt[t],GRB_EQUAL, nt, constraint_name);
	for (int k=0; k<dt.K; k++){
	  constraint_name = "set_var_Nkt_k="+to_string(k)+"t="+to_string(t);
	  md.addConstr(var.Nkt[t*dt.K+k],GRB_EQUAL, nkt[k], constraint_name);
	}
      }
    }
  }
}

int build_model::compute_number_branchings(){
  int nb = 0;
  if (mt.univ){
    for (int t=0; t<param.N; t++){
      for (int j= 0; j<param.J; j++){
	nb += (int)(var.a[t*param.J+j].get(GRB_DoubleAttr_X) >=0.99); 
      }
    }
  }
  else{
    for (int t=0; t<param.N; t++){
      for (int j= 0; j<param.J; j++){
	nb += (int)(var.s[t*param.J+j].get(GRB_DoubleAttr_X) >=0.99); 
      }
    }
  }
  return nb;
}

void build_model::buildTree(Tree& T){ 
  for (int t = 0; t < param.N; t++){
    T.b[t] = var.b[t].get(GRB_DoubleAttr_X);
    if (mt.eps){
      T.eps[t] = var.eps[t].get(GRB_DoubleAttr_X);
    }
    else{
      T.eps[t] = 0.0;
    }
    for (int j=0; j < param.J; j++){
      T.a[t*param.J+j] = var.a[t*param.J+j].get(GRB_DoubleAttr_X);
    }
  }
  
  if (mt.base == baseModel::F){
    for (int t=0; t<param.N+param.L; t++){ // initialisation, on met les feuilles à -2
      T.c[t] = -1 + -1*(t>=param.N);
    }
    
    for (int t=0; t<param.N; t++){
      int clas = -1;
      for (int k=0; k<param.K; k++){
	if (var.c[t*param.K+k].get(GRB_DoubleAttr_X) >= 0.99){
	  clas = k;
	  break;
	}
      }
      if (clas != -1 && T.c[T.right_most_leaf[t]] == -2){
	if (t<param.N){
	  for (auto l: T.leaves[t]){
	    T.c[l] = -1;
	  }
	}
	T.c[t] = clas;
      }
    }

    for (int t=param.N; t<param.N+param.L; t++){
      int clas = -1;
      for (int k=0; k<param.K; k++){
	if (var.c[t*param.K+k].get(GRB_DoubleAttr_X) >= 0.99){
	  clas = k;
	  break;
	}
      }
      if (clas != -1 && T.c[t] == -2){
	T.c[t] = clas;
      }
    }
  }
  else{ // leaf_class tables are the same but:
    int* leaf_class_1 = new int[param.L]; // it wont change in the process
    int* leaf_class_2 = new int[param.L]; // it will change in the process
    for (int t = 0; t < param.L; t++){
      leaf_class_1[t] = -1;
      leaf_class_2[t] = -1;
      
      int A_t = get_A(t,param.N);
      bool is_real_leaf = true;
      if (t != param.L - 1){
	is_real_leaf = var.d[A_t].get(GRB_DoubleAttr_X) >= 0.99;
      }
      for (int k=0; k<param.K; k++){
	
	if (var.c[t*param.K+k].get(GRB_DoubleAttr_X) >= 0.99 and is_real_leaf){
	  leaf_class_1[t] = k;
	  leaf_class_2[t] = k;
	  break;
	}
      }
    }
    // cas particulier du root node
    if (leaf_class_1[T.right_most_leaf[1]-param.N] >= 0 and var.d[0].get(GRB_DoubleAttr_X)==0.0){
      T.c[0] = leaf_class_1[T.right_most_leaf[0]-param.N];
      leaf_class_2[T.right_most_leaf[0]-param.N] = -1;
    }
    else{
      T.c[0] = -1;
    }
    for (int t = 1; t<param.N; t++){
      int a_t = (t + 1)/2-1,
	l_t = (t+1)*2-1;
      int l_a_t = (a_t+1)*2-1,
	r_a_t = (a_t+1)*2;
      bool ancestor_branching = (leaf_class_1[T.right_most_leaf[l_a_t]-param.N]>=0)&&(leaf_class_1[T.right_most_leaf[r_a_t]-param.N]>=0),
	current_node_potential_class = leaf_class_1[T.right_most_leaf[l_t]-param.N] == -1;
      if (ancestor_branching&&current_node_potential_class){ // if the ancestor of t is branching and the right most leaf of t's left son does not assign a class, it means t assigns a class
	T.c[t] = leaf_class_1[T.right_most_leaf[t]-param.N]; // the class t assigns is the class of t's right most leaf
	leaf_class_2[T.right_most_leaf[t]-param.N] = -1; // in the second table, we specifiy that t's right most leaf is not actually assigning a class
      }
      else{
	T.c[t] = -1;
      }
    }
    for (int t = 0; t<param.L; t++){ // if a leaf was a "proxy" leaf for an ancestor then before we would have put -1 otherwise it actually assign the class given by OCT, therefore we can just copy leaf_class_2
      T.c[t+param.N] = leaf_class_2[t];
    }
  }
  /*
  if (!mt.univ){
    T.correctSplitFunctions();
    }*/
}


void build_model::get_z(int z[]){
  
  for (int i = 0; i<param.I;i++){
    for (int t = 0; t<param.L; t++){
      z[i*param.L+t] = round(var.z[i*param.L+t].get(GRB_DoubleAttr_X));
    }
  }

  
  /*
  for (int k = 0; k<param.K;k++){
    for (int t = 0; t<param.N+param.L; t++){
      z[t*param.K+k] = var.c[t*param.K+k].get(GRB_DoubleAttr_X);
    }
  }
  */
  /*
  for (int t = 0; t<param.L; t++){
    z[t] = var.l[t].get(GRB_DoubleAttr_X);
  }
  */

  /*
  for (int j = 0; j<param.J;j++){
    for (int t = 0; t<param.N; t++){
      z[t*param.J+j] = var.a[t*param.J+j].get(GRB_DoubleAttr_X);
    }
  }
  /*
  
  /*
  for (int i=0; i<param.I; i++){
    for (int e=0; e<param.E; e++){
      z[i*param.E+e] = var.u[i*param.E+e].get(GRB_DoubleAttr_X);
    }
  }
  */
}

solution build_model::solve(GRBModel& md, double time_limit){
  md.set("TimeLimit", to_string(time_limit));
  
  freopen("gurobi_text.txt", "a", stdout);
  md.optimize();
  freopen("/dev/tty", "w", stdout);
  remove("gurobi_text.txt");

  double obj_val = md.get(GRB_DoubleAttr_ObjVal);
  int nb_br = compute_number_branchings();

  int err_tr = 0;
  err_tr = ceil(obj_val - param.alph*nb_br)*param.L_hat;
  
  double rr = read_root_rel_objvalue();
  
  //int* z = new int[param.J*param.N];
  //get_z(z);
  
  Tree T = Tree(param.D,param.J,param.K);
  buildTree(T);

  //double gap = 0.0;
  double gap = md.get(GRB_DoubleAttr_MIPGap); // ya un bug "Unable to retrieve attribute 'MIPGap'" donc j'enlève momentanément

  return solution(T,obj_val,err_tr,nb_br,md.get(GRB_DoubleAttr_Runtime),md.get(GRB_DoubleAttr_MIPGap),(int)md.get(GRB_DoubleAttr_NodeCount),rr);   
}

double build_model::read_root_rel_objvalue(){
  /*
  fstream file;
  file.open("gurobi_text.txt",ios::in);
  string line;
  file>>line;
  while (line != "Root"){
    file >> line;
  }
  file >> line;
  file >> line;
  double obj;
  file >> obj;
  file.close();
  remove("gurobi_text.txt");
  */

  //return obj;
  
  return 0;
}
