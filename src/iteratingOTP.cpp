#include "iteratingOTP.h"

void solClust::write(string filename){
  fstream file;
  //file.open(filename,ios::out | std::ofstream::trunc);
  file.open(filename,ios::out);
  file << totTime << endl;
  file << clTime << endl;
  file << nbClInitial <<endl;
  file << nbClFinal << endl;
  file << errTr << endl;
  file << errTst << endl;
  file << nbIter << endl;

  for (int i=0; i<nbIter; i++){
    file << solveTime[i] << endl;
    file << centeringTime[i] << endl;
    file << isOpti[i] << endl;
    file << nbCl[i] << endl;
    file << errCl[i] << endl;
    file << errDt[i] << endl;
  }

  file.close();
}

bool isIntersecting(Tree T, clustering& cl, dataset& initialDt, dataset& currentDt, bool univ, bool updateCl){
  vector<int>* IL = new vector<int>[T.N];
  vector<int>* IR = new vector<int>[T.N];
  T.compute_ILIR(IL, IR, currentDt, true);
  
  vector<int>* isCut = new vector<int>[cl.clusters.size()];

  bool B = false;
  GRBEnv env = GRBEnv();
  for (int t=0; t<T.N; t++){    
    if (IL[t].size() !=0 && IR[t].size() !=0){
      int sum_sj_star = 0;
      for (int j=0; j<initialDt.J; j++){
	if (!equalTo(T.a[t*initialDt.J+j],0)){
	  sum_sj_star += 1;
	}
      }

      int nbPoints = 0;
      for (int i=0; i<IL[t].size(); i++){
	nbPoints += cl.clusters[IL[t][i]].pts.size();
      }
      for (int i=0; i<IR[t].size(); i++){
	nbPoints += cl.clusters[IR[t][i]].pts.size();
      }

      GRBModel m = GRBModel(env);

      GRBVar b1_var;
      b1_var = m.addVar(-1.0, 1.0, 0.0, GRB_CONTINUOUS, "b1");
      GRBVar b2_var;
      b2_var = m.addVar(-1.0, 1.0, 0.0, GRB_CONTINUOUS, "b2");
      m.addConstr(b1_var+mu,GRB_LESS_EQUAL,b2_var,"pour_eviter_la_triche");
      
      GRBVar* a_var;
      double* lb_a = new double[initialDt.J];
      double* ub_a = new double[initialDt.J];
      double* coef_a = new double[initialDt.J];
      char* type_a = new char[initialDt.J];
      string* name_a = new string[initialDt.J];
      for (int j=0; j<initialDt.J; j++){
	lb_a[j] = -1.0;
	ub_a[j] = 1.0;
	coef_a[j] = 0.0;
	type_a[j] = GRB_CONTINUOUS;
	name_a[j] = "a" + to_string(j);
      }
      a_var = m.addVars(lb_a,ub_a,coef_a,type_a,name_a,initialDt.J);

      GRBVar* ah_var;
      double* lb_ah = new double[initialDt.J];
      double* ub_ah = new double[initialDt.J];
      double* coef_ah = new double[initialDt.J];
      char* type_ah = new char[initialDt.J];
      string* name_ah = new string[initialDt.J];
      for (int j=0; j<initialDt.J; j++){
	lb_ah[j] = 0.0;
	ub_ah[j] = 1.0;
	coef_ah[j] = 0.0;
	type_ah[j] = GRB_CONTINUOUS;
	name_ah[j] = "ah" + to_string(j);
      }
      ah_var = m.addVars(lb_ah,ub_ah,coef_ah,type_ah,name_ah,initialDt.J);

      GRBVar* s_var;
      s_var = m.addVars(initialDt.J, GRB_BINARY);


      GRBVar* psi_var;
      psi_var = m.addVars(currentDt.I, GRB_BINARY);

      GRBVar* z_var;
      z_var = m.addVars(nbPoints*2, GRB_BINARY);

      string constraint_name;
      GRBLinExpr sum_sj = 0;
      GRBLinExpr sum_ahj = 0;
      for (int j=0; j<initialDt.J;j++){ // -sj <= aj <= sj et -ah_j <= aj <= ah_j
	constraint_name = "a_j_geq_-sj_j="+to_string(j);
	m.addConstr(a_var[j],GRB_GREATER_EQUAL,-s_var[j],constraint_name);
	constraint_name = "aj_leq_sj_j="+to_string(j);
	m.addConstr(a_var[j],GRB_LESS_EQUAL,s_var[j],constraint_name);
	sum_sj += s_var[j];
	
	constraint_name = "a_j_geq_-ahj_j="+to_string(j);
	m.addConstr(a_var[j],GRB_GREATER_EQUAL,-ah_var[j],constraint_name);
	constraint_name = "aj_leq_ahj_j="+to_string(j);
	m.addConstr(a_var[j],GRB_LESS_EQUAL,ah_var[j],constraint_name);
	sum_ahj += ah_var[j];
      }
      constraint_name = "lesser_or_equal_complexity"; // sum sj <= sum s*j
      m.addConstr(sum_sj,GRB_LESS_EQUAL,sum_sj_star,constraint_name);

      constraint_name = "limit_sum_ahj"; // sum ahj <= 1
      m.addConstr(sum_ahj,GRB_LESS_EQUAL,1,constraint_name);

      for (int i=0; i<nbPoints; i++){ // zi,l + zi,r = 1
	constraint_name = "left_or_right_"+to_string(i);
	m.addConstr(z_var[i*2]+z_var[i*2+1],GRB_EQUAL,1,constraint_name);	
      }

      int cpt=0;
      for (int g=0; g<IL[t].size(); g++){
	GRBLinExpr sum_ax = 0;
	for (int j=0; j<initialDt.J; j++){
	  sum_ax += a_var[j]*currentDt.X[IL[t][g]*currentDt.J+j];
	}
	constraint_name = "representative_g_goes_into_left_branch_g="+to_string(g);
	m.addConstr(sum_ax,GRB_LESS_EQUAL,b1_var,constraint_name);

	for (int i=0; i<cl.clusters[IL[t][g]].pts.size(); i++){
	  GRBLinExpr sum_ax_i = 0;
	  for (int j=0; j<initialDt.J; j++){
	    sum_ax_i += a_var[j]*initialDt.X[cl.clusters[IL[t][g]].pts[i]*initialDt.J+j];
	  }
	  constraint_name = "going_to_left_branch_i="+to_string(cpt+i);
	  m.addConstr(sum_ax_i,GRB_LESS_EQUAL,b1_var+2*(1-z_var[(cpt+i)*2]),constraint_name);
	  constraint_name = "going_to_right_branch_i="+to_string(cpt+i);
	  m.addConstr(sum_ax_i,GRB_GREATER_EQUAL,b2_var-2*(1-z_var[(cpt+i)*2+1]),constraint_name);

	  constraint_name = "cluster_cassé?_g="+to_string(g)+"_i"+to_string(i);
	  m.addConstr(psi_var[g],GRB_GREATER_EQUAL, z_var[(cpt+i)*2+1], constraint_name);
	}	
	cpt += cl.clusters[IL[t][g]].pts.size();
      }
      for (int g=0; g<IR[t].size(); g++){
	GRBLinExpr sum_ax = 0;
	for (int j=0; j<initialDt.J; j++){
	  sum_ax += a_var[j]*currentDt.X[IR[t][g]*currentDt.J+j];
	}
	constraint_name = "representative_g_goes_into_right_branch_g="+to_string(g+IL[t].size());
	m.addConstr(sum_ax,GRB_GREATER_EQUAL,b2_var,constraint_name);
	
	for (int i=0; i<cl.clusters[IR[t][g]].pts.size(); i++){
	  GRBLinExpr sum_ax_i = 0;
	  for (int j=0; j<initialDt.J; j++){
	    sum_ax_i += a_var[j]*initialDt.X[cl.clusters[IR[t][g]].pts[i]*initialDt.J+j];
	  }
	  
	  constraint_name = "going_to_left_branch_i="+to_string(cpt+i);
	  m.addConstr(sum_ax_i,GRB_LESS_EQUAL,b1_var+2*(1-z_var[(cpt+i)*2]),constraint_name);
	  constraint_name = "going_to_right_branch_i="+to_string(cpt+i);
	  m.addConstr(sum_ax_i,GRB_GREATER_EQUAL,b2_var-2*(1-z_var[(cpt+i)*2+1]),constraint_name);

	  constraint_name = "cluster_cassé?_g="+to_string(g+IL[t].size())+"_i"+to_string(i);
	  m.addConstr(psi_var[g+IL[t].size()],GRB_GREATER_EQUAL, z_var[(cpt+i)*2], constraint_name);
	}
	
	cpt += cl.clusters[IR[t][g]].pts.size();
      }
      GRBLinExpr objFunc = 0;
      for (int g=0; g<IL[t].size()+IR[t].size(); g++){
	objFunc += psi_var[g]; // on minimise le nomber de clusters cassés
      }
      objFunc += (b1_var - b2_var)/4; // et la marge faite à chaque séparation
      m.setObjective(objFunc,GRB_MINIMIZE);

      if (setPrecision){
	m.set(GRB_DoubleParam_IntFeasTol, GLOBAL_PRECISION); // val par def : 10-5, on passe à 10-6
	m.set(GRB_DoubleParam_FeasibilityTol, GLOBAL_PRECISION); // comme la valeur par defaut
	m.set(GRB_DoubleParam_OptimalityTol, GLOBAL_PRECISION); // comme la valeur par defaut
      }

      //m.write("model" + to_string(t)+".lp");
      freopen("gurobi_text.txt", "a", stdout);

      m.optimize();
      
      freopen("/dev/tty", "w", stdout);
      // updating the tree

      T.b[t] = (b1_var.get(GRB_DoubleAttr_X) + b2_var.get(GRB_DoubleAttr_X))/ 2;
      for (int j=0; j<initialDt.J; j++){
	T.a[t*initialDt.J+j] = a_var[j].get(GRB_DoubleAttr_X);
      }

      if (greaterThan(m.get(GRB_DoubleAttr_ObjVal),0)){
	B = true;
	for (int g=0; g<IL[t].size(); g++){
	  if (psi_var[g].get(GRB_DoubleAttr_X) > 0){
	    isCut[IL[t][g]].push_back(t);
	  }
	}
	for (int g=0; g<IR[t].size(); g++){
	  if (psi_var[g+IL[t].size()].get(GRB_DoubleAttr_X) > 0){
	    isCut[IR[t][g]].push_back(t);
	  }
	}
      }      
    }
  }

  if (updateCl and B){
    int* leavesOfAllPoints = new int[initialDt.I];
    T.predict_leaves(initialDt, leavesOfAllPoints);
    for (int g=0; g<currentDt.I; g++){
      if (isCut[g].size()>0){
	vector<vector<int>> division;
	map<int,int> leafOfDiv;
	for (int i=0; i<cl.clusters[g].pts.size(); i++){
	  if (leafOfDiv.count(leavesOfAllPoints[cl.clusters[g].pts[i]]) > 0){
	    division[leafOfDiv[leavesOfAllPoints[cl.clusters[g].pts[i]]]].push_back(cl.clusters[g].pts[i]);
	  }
	  else{
	    leafOfDiv[leavesOfAllPoints[cl.clusters[g].pts[i]]] = division.size();
	    division.push_back({cl.clusters[g].pts[i]});
	  }
	}
	cl.breakCluster(g,division,initialDt);
      }
    }
  }

  return B;
}

solClust iteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, int timeL){
  clustering clust = cl.clustering_copy();
  parameters pm = p.parameters_copy();

  bool keepGoing = true;

  solClust solC = solClust();
  solC.nbClInitial = clust.clusters.size();

  //clust.showClusters();
  
  GRBEnv env = GRBEnv();
  while(keepGoing and clust.clusters.size() < dt.I){
    time_t t1 = time (NULL);
    
    dataset currentDt = clust.createDt(dt);
    pm.update(currentDt);
    GRBModel md = GRBModel(env);
    build_model model = build_model(md,currentDt,modelt,pm);
    solution sol = model.solve(md,timeL-solC.totTime);
    
    time_t t2 = time (NULL);

    keepGoing = isIntersecting(sol.T, clust, dt, currentDt, modelt.univ, true);
    time_t t3 = time (NULL);

    solC.T = sol.T;
    //sol.T.write_tree("arbre"+to_string(solC.nbIter)+".txt");
    solC.solveTime.push_back(t2-t1);
    solC.centeringTime.push_back(t3-t2);
    solC.totTime += t3 - t1;
    solC.nbIter += 1;
    solC.isOpti.push_back(sol.time < timeL-solC.totTime);
    solC.nbCl.push_back(clust.clusters.size());
    solC.errCl.push_back(sol.T.prediction_errors(currentDt));
    solC.errDt.push_back(sol.T.prediction_errors(dt));

    if (solC.totTime > timeL){
      break;
    }
  }
  solC.nbClFinal = clust.clusters.size();

  solC.finalCl = clust.clustering_copy();

  //clust.showClusters();

  return solC;
}
