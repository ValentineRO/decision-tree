#include "iteratingOTP.h"

void solClust::write(string filename){
  fstream file;
  //file.open(filename,ios::out | std::ofstream::trunc);
  file.open(filename,ios::out);
  file << totTime << endl;
  file << clTime << endl;
  file << nbClInitial <<endl;
  file << nbClFinal << endl;
  file << pseudoGap << endl;
  file << nbIter << endl;

  for (int i=0; i<nbIter; i++){
    file << solveTime[i] << endl;
    file << centeringTime[i] << endl;
    file << obj[i] << endl;
    file << nbCl[i] << endl;
    file << errCl[i] << endl;
    file << errDt[i] << endl;
  }

  file.close();
}

bool isIntersecting(Tree T, clustering& cl, dataset& initialDt, dataset& currentDt, bool univ, int Nmin){
  vector<int>* IL = new vector<int>[T.N];
  vector<int>* IR = new vector<int>[T.N];
  T.compute_ILIR(IL, IR, currentDt, true);
  
  bool B = false;
  GRBEnv env = GRBEnv();
  for (int t=0; t<T.N; t++){
    if (T.c[t] == -1 and IL[t].size()!=0 and IR[t].size()!=0){
      // calcul du nombre de points au noeud t (c'est oké car IL, IR et cl sont updatés)
      int nbPoints = 0;
      for (int i=0; i<IL[t].size(); i++){
	//cout << "a gauche " << IL[t][i] << endl;
	nbPoints += cl.clusters[IL[t][i]].pts.size();
      }
      for (int i=0; i<IR[t].size(); i++){
	//cout << "a droite " << IR[t][i] << endl;
	nbPoints += cl.clusters[IR[t][i]].pts.size();
      }

      GRBModel m = GRBModel(env);

      GRBVar* psi_var;
      psi_var = m.addVars(IL[t].size()+IR[t].size(), GRB_BINARY);
      GRBVar* z_var;
      z_var = m.addVars(nbPoints*2, GRB_BINARY);
      GRBVar ll_var;
      ll_var = m.addVar(0.0, 1.0, 0.0, GRB_BINARY, "ll");
      GRBVar lr_var;
      lr_var = m.addVar(0.0, 1.0, 0.0, GRB_BINARY, "lr");

      GRBVar b1_var, b2_var;
      GRBVar* a_var;
      GRBVar* ah_var;
      GRBVar* s_var;
      if (univ){
	b1_var = m.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "b1");
	b2_var = m.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "b2");
	m.addConstr(b1_var,GRB_LESS_EQUAL,b2_var-mu*ll_var,"b1<=b2-d*mu");

	double* lb_a = new double[initialDt.J];
	double* ub_a = new double[initialDt.J];
	double* coef_a = new double[initialDt.J];
	char* type_a = new char[initialDt.J];
	string* name_a = new string[initialDt.J];
	for (int j=0; j<initialDt.J; j++){
	  lb_a[j] = 0.0;
	  ub_a[j] = 1.0;
	  coef_a[j] = 0.0;
	  type_a[j] = GRB_BINARY;
	  name_a[j] = "a" + to_string(j);
	}
	a_var = m.addVars(lb_a,ub_a,coef_a,type_a,name_a,initialDt.J);

	GRBLinExpr sum_aj = 0;
	for (int j=0; j<T.J; j++){
	  sum_aj += a_var[j];
	}
	m.addConstr(sum_aj,GRB_EQUAL,ll_var,"sum_a=d");
      }
      else{
	 b1_var = m.addVar(-1.0, 1.0, 0.0, GRB_CONTINUOUS, "b1");
	 b2_var = m.addVar(-1.0, 1.0, 0.0, GRB_CONTINUOUS, "b2");
	 m.addConstr(b1_var,GRB_LESS_EQUAL,b2_var-mu*ll_var,"b1<=b2");

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

	 s_var = m.addVars(initialDt.J, GRB_BINARY);

	 GRBLinExpr sum_sj, sum_ahj;

	 for (int j=0; j<T.J; j++){
	   m.addConstr(-s_var[j],GRB_LESS_EQUAL,a_var[j],"-sj<=aj_j="+to_string(j));
	   m.addConstr(s_var[j],GRB_GREATER_EQUAL,a_var[j],"-sj>=aj_j="+to_string(j));
	   m.addConstr(-ah_var[j],GRB_LESS_EQUAL,a_var[j],"-ah_j<=aj_j="+to_string(j));
	   m.addConstr(ah_var[j],GRB_GREATER_EQUAL,a_var[j],"-ah_j>=aj_j="+to_string(j));

	   sum_sj += s_var[j];
	   sum_ahj += ah_var[j];
	 }
	 m.addConstr(sum_ahj,GRB_LESS_EQUAL,ll_var,"sum_ahj<=d");
	 
	 int sum_sj_star = 0;
	 for (int j=0; j<initialDt.J; j++){
	   if (!equalTo(T.a[t*initialDt.J+j],0)){
	     sum_sj_star += 1;
	   }
	 }
	 m.addConstr(sum_sj,GRB_LESS_EQUAL,sum_sj_star,"sum_sj<=s_star");
	 m.addConstr(sum_sj,GRB_GREATER_EQUAL,ll_var,"sum_sj>=d");
	 m.addConstr(b1_var,GRB_GREATER_EQUAL,-ll_var,"-d<=b1");
      }
      m.addConstr(b2_var,GRB_LESS_EQUAL,ll_var,"b2<=d");

      vector<int> cumulL, cumulR;
      int cpt=0;

      string constraint_name;
      for (int g=0; g<IL[t].size(); g++){
	GRBLinExpr sum_ax = 0;
	for (int j=0; j<initialDt.J; j++){
	  sum_ax += a_var[j]*cl.clusters[IL[t][g]].bar[j];
	}
        
	constraint_name = "representative_g_goes_into_left_branch_g="+to_string(g);
	m.addConstr(sum_ax,GRB_LESS_EQUAL,b1_var,constraint_name);

	for (int i=0; i<cl.clusters[IL[t][g]].pts.size(); i++){
	  GRBLinExpr sum_ax_i = 0;
	  for (int j=0; j<initialDt.J; j++){
	    sum_ax_i += a_var[j]*initialDt.X[cl.clusters[IL[t][g]].pts[i]*initialDt.J+j];
	  }

	  if (univ){
	    constraint_name = "going_to_left_branch_i="+to_string(cpt+i);
	    m.addConstr(sum_ax_i,GRB_LESS_EQUAL,b1_var+(1+mu)*(1-z_var[(cpt+i)*2]),constraint_name);
	    constraint_name = "going_to_right_branch_i="+to_string(cpt+i);
	    m.addConstr(sum_ax_i,GRB_GREATER_EQUAL,b2_var-(1-z_var[(cpt+i)*2+1]),constraint_name);
	  }
	  else{
	    constraint_name = "going_to_left_branch_i="+to_string(cpt+i);
	    m.addConstr(sum_ax_i,GRB_LESS_EQUAL,b1_var+(2+mu)*(1-z_var[(cpt+i)*2]),constraint_name);
	    constraint_name = "going_to_right_branch_i="+to_string(cpt+i);
	    m.addConstr(sum_ax_i,GRB_GREATER_EQUAL,b2_var-2*(1-z_var[(cpt+i)*2+1]),constraint_name);
	  }

	  constraint_name = "cluster_cassé?_g="+to_string(g)+"_i"+to_string(i);
	  m.addConstr(psi_var[g],GRB_GREATER_EQUAL, z_var[(cpt+i)*2+1], constraint_name);
	}
	cumulL.push_back(cpt);
	cpt += cl.clusters[IL[t][g]].pts.size();
      }
      
      for (int g=0; g<IR[t].size(); g++){
	GRBLinExpr sum_ax = 0;
	for (int j=0; j<initialDt.J; j++){
	  sum_ax += a_var[j]*cl.clusters[IR[t][g]].bar[j];
	}
	constraint_name = "representative_g_goes_into_right_branch_g="+to_string(g+IL[t].size());
	m.addConstr(sum_ax,GRB_GREATER_EQUAL,b2_var,constraint_name);
	
	for (int i=0; i<cl.clusters[IR[t][g]].pts.size(); i++){
	  GRBLinExpr sum_ax_i = 0;
	  for (int j=0; j<initialDt.J; j++){
	    sum_ax_i += a_var[j]*initialDt.X[cl.clusters[IR[t][g]].pts[i]*initialDt.J+j];
	  }

	  if (univ){
	    constraint_name = "going_to_left_branch_i="+to_string(cpt+i);
	    m.addConstr(sum_ax_i,GRB_LESS_EQUAL,b1_var+(1+mu)*(1-z_var[(cpt+i)*2]),constraint_name);
	    constraint_name = "going_to_right_branch_i="+to_string(cpt+i);
	    m.addConstr(sum_ax_i,GRB_GREATER_EQUAL,b2_var-(1-z_var[(cpt+i)*2+1]),constraint_name);
	  }
	  else{
	    constraint_name = "going_to_left_branch_i="+to_string(cpt+i);
	    m.addConstr(sum_ax_i+mu,GRB_LESS_EQUAL,b1_var+(2+mu)*(1-z_var[(cpt+i)*2]),constraint_name);
	    constraint_name = "going_to_right_branch_i="+to_string(cpt+i);
	    m.addConstr(sum_ax_i,GRB_GREATER_EQUAL,b2_var-2*(1-z_var[(cpt+i)*2+1]),constraint_name);
	  }
	  
	  constraint_name = "cluster_cassé?_g="+to_string(g+IL[t].size())+"_i"+to_string(i);
	  m.addConstr(psi_var[g+IL[t].size()],GRB_GREATER_EQUAL, z_var[(cpt+i)*2], constraint_name);
	}
	cumulR.push_back(cpt);
	cpt += cl.clusters[IR[t][g]].pts.size();
      }

      GRBLinExpr sum_zil = 0,
	sum_zir = 0;
      for (int i=0; i<nbPoints; i++){
	sum_zil += z_var[i*2];
	sum_zir += z_var[i*2+1];
	m.addConstr(z_var[i*2],GRB_LESS_EQUAL,ll_var,"zil<=ll_i="+to_string(i));
	m.addConstr(z_var[i*2+1],GRB_LESS_EQUAL,lr_var,"zir<=lr_i="+to_string(i));
	m.addConstr(z_var[i*2] + z_var[i*2+1],GRB_EQUAL,1,"zil+zir=1"+to_string(i));
      }
      m.addConstr(sum_zil,GRB_GREATER_EQUAL,Nmin*ll_var,"sum_zil>=Nmin*ll");
      m.addConstr(sum_zir,GRB_GREATER_EQUAL,Nmin*lr_var,"sum_zir>=Nmin*lr");
      
      GRBLinExpr objFunc = 0;
      for (int g=0; g<IL[t].size()+IR[t].size(); g++){
	objFunc += psi_var[g]; // on minimise le nomber de clusters cassés
      }
      objFunc += (b1_var - b2_var)/4; // et la marge faite à chaque séparation
      objFunc += (IL[t].size()+IR[t].size())*(1-ll_var); // on pénalise beaucoup le fait de devenir une feuille
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
      remove("gurobi_text.txt");

      if (equalTo(ll_var.get(GRB_DoubleAttr_X),0)){ // on transforme le noeud en feuille
	int* nbPerClass = new int[T.K];
	for (int k=0; k<T.K; k++){
	  nbPerClass[k] = 0;
	}
	for (int c=0; c<IL[t].size(); c++){
	  nbPerClass[cl.clusters[IL[t][c]].lb] += cl.clusters[IL[t][c]].pts.size();
	}
	for (int c=0; c<IR[t].size(); c++){
	  nbPerClass[cl.clusters[IR[t][c]].lb] += cl.clusters[IR[t][c]].pts.size();
	}
	T.c[t] = 0;
	for (int k=1; k<T.K; k++){
	  if(nbPerClass[k] > nbPerClass[T.c[t]]){
	    T.c[t] = k;
	  }
	}
	T.b[t] = 0;
	for (int j=0; j<T.J; j++){
	  T.a[t*T.J+j] = 0;
	}
	vector<int> changingNodes = {t*2+1,t*2+2};
	while(changingNodes.size() > 0){
	  int currentNode = changingNodes.back();
	  changingNodes.pop_back();
	  T.c[currentNode] == -1;
	  if (currentNode < T.N){
	    for (int j=0; j<T.J; j++){
	      T.a[currentNode*T.J+j] = 0;
	    }
	    T.b[currentNode] = 0;
	    IL[currentNode] = {};
	    IR[currentNode] = {};
	    changingNodes.push_back(currentNode*2+1);
	    changingNodes.push_back(currentNode*2+2);
	  }	    
	}
      }
      else{
	T.b[t] = (b1_var.get(GRB_DoubleAttr_X) + b2_var.get(GRB_DoubleAttr_X))/ 2;
	for (int j=0; j<initialDt.J; j++){
	  T.a[t*initialDt.J+j] = a_var[j].get(GRB_DoubleAttr_X);
	}

	int currentILsize = IL[t].size(),
	  currentIRsize = IR[t].size();
	
	if (greaterThan(m.get(GRB_DoubleAttr_ObjVal),0)){
	  B = true;
	  for (int g=0; g<currentILsize; g++){
	    if (psi_var[g].get(GRB_DoubleAttr_X) > 0){
	      vector<vector<int>> division;
	      division.push_back({});
	      division.push_back({});
	      for (int i=0; i<cl.clusters[IL[t][g]].pts.size(); i++){
		if (z_var[(cumulL[g]+i)*2+1].get(GRB_DoubleAttr_X) > 0){
		  division[1].push_back(cl.clusters[IL[t][g]].pts[i]);
		}
		else{
		  division[0].push_back(cl.clusters[IL[t][g]].pts[i]);
		}
	      }
	      cl.breakCluster(IL[t][g], division, initialDt);
	      int lf = T.predict_leaf(cl.clusters.back().bar);
	      while(lf>t){
		int alf = (lf-1)/2;
		if(lf<T.N){
		  if (lf== alf*2+1){
		    IL[alf].push_back(cl.clusters.size()-1);
		  }
		  else{
		    IR[alf].push_back(cl.clusters.size()-1);
		  }
		}
		lf = alf;
	      }
	    }
	  }
	  for (int g=0; g<currentIRsize; g++){
	    if (psi_var[g+currentILsize].get(GRB_DoubleAttr_X) > 0){
	      vector<vector<int>> division;
	      division.push_back({});
	      division.push_back({});
	      for (int i=0; i<cl.clusters[IR[t][g]].pts.size(); i++){
		if (z_var[(cumulR[g]+i)*2].get(GRB_DoubleAttr_X) > 0){
		  division[1].push_back(cl.clusters[IR[t][g]].pts[i]);
		}
		else{
		  division[0].push_back(cl.clusters[IR[t][g]].pts[i]);
		}
	      }
	      cl.breakCluster(IR[t][g], division, initialDt);
	      int lf = T.predict_leaf(cl.clusters.back().bar);
	      while(lf>t){
		int alf = (lf-1)/2;
		if(lf<T.N){
		  if (lf== alf*2+1){
		    IL[alf].push_back(cl.clusters.size()-1);
		  }
		  else{
		    IR[alf].push_back(cl.clusters.size()-1);
		  }
		}
		lf = alf;
	      }
	    }
	  }
	}
      }
    }
  }
  
  return B;
}
/*
bool isIntersectingOG(Tree T, clustering& cl, dataset& initialDt, dataset& currentDt, bool univ, bool updateCl, int Nmin){
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

      GRBLinExpr sum_left = 0,
	sum_right = 0;
      for (int i=0; i<nbPoints; i++){
	sum_left += z_var[i*2];
	sum_right += z_var[i*2+1];
      }
      constraint_name = "verif_Nmin_left";
      m.addConstr(sum_left, GRB_GREATER_EQUAL, Nmin, constraint_name);
      constraint_name = "verif_Nmin_right";
      m.addConstr(sum_right, GRB_GREATER_EQUAL, Nmin, constraint_name);
      
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
      remove("gurobi_text.txt");

      
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
*/

solClust approxIteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, double timeL){
  clustering clust = cl.clustering_copy();
  parameters pm = p.parameters_copy();

  pm.MIPFocus = 1;

  bool keepGoing = true;

  solClust solC = solClust();
  solC.nbClInitial = clust.clusters.size();

  //clust.showClusters();
  int cpt = 0;
  
  GRBEnv env = GRBEnv();
  while(timeL-solC.totTime> 10.0 and keepGoing and clust.clusters.size() < dt.I){
    cpt += 1;

    time_t t1 = time (NULL);
    
    dataset currentDt = clust.createDt(dt);
    pm.update(currentDt);

    GRBModel md = GRBModel(env);
    build_model model = build_model(md,currentDt,modelt,pm);

    Tree cartTree = pythonCART(currentDt,pm.D,pm.C,pm.Nmin);
    if (cpt == 1){
      model.add_warmstart(cartTree, currentDt);
    }
    else{
      if (cartTree.prediction_errors(currentDt) > solC.T.back().prediction_errors(currentDt)){
	model.add_warmstart(solC.T.back(), currentDt);
      }
      else{
	model.add_warmstart(cartTree, currentDt);
      }
    }

    solution sol = model.solve2(md,timeL-solC.totTime);
    time_t t2 = time (NULL);
    solC.nbCl.push_back(clust.clusters.size()); // on met le nombre de cl considérés avant de potentiellement changer les clusters
    keepGoing = isIntersecting(sol.T, clust, dt, currentDt, modelt.univ, pm.Nmin);
    time_t t3 = time (NULL);

    /*
    // on récupère des infos mais c'est à court terme
    freopen("bsMIPFocus1.txt", "a", stdout);
    cout << "Iteration " << cpt << endl;
    cout << cartTree.prediction_errors(currentDt) << endl;
    cout << sol.obj << endl;
    cout << sol.gap << endl;
    cout << sol.T.prediction_errors(dt) << endl;
    cout << sol.time << endl;
    freopen("/dev/tty", "w", stdout);
    // a delete plus tard */

    solC.T.push_back(sol.T);
    solC.solveTime.push_back(t2-t1);
    solC.centeringTime.push_back(t3-t2);
    solC.gap.push_back(sol.gap);
    solC.totTime += t3 - t1;
    solC.nbIter += 1;
    solC.obj.push_back(sol.obj);
    solC.errCl.push_back(solC.T.back().prediction_errors(currentDt));
    solC.errDt.push_back(solC.T.back().prediction_errors(dt));

    if (solC.totTime > timeL){
      break;
    }
  }
  solC.nbClFinal = clust.clusters.size();
  solC.finalCl = clust.clustering_copy();

  if (solC.errCl.back() == 0){
    if (solC.errDt.back() == 0){
      solC.pseudoGap = 0.0;
    }
    else{
      solC.pseudoGap = 1.0;
    }
  }
  else{
    solC.pseudoGap = (double)(solC.errDt.back()-solC.errCl.back())/(double)solC.errDt.back();
  }
  
  return solC;
}

solClust iteratingCART(dataset& dt, clustering& cl, parameters p){
  clustering clust = cl.clustering_copy();

  bool keepGoing = true;

  solClust solC = solClust();
  solC.nbClInitial = clust.clusters.size();

  Tree OGcartTree = pythonCART(dt, p.D,p.C,p.Nmin);

  //OGcartTree.write_tree("iteration0.txt");

  cout << "Nombre de données initial " << dt.I << endl;
  cout << "Valeur CART initiale " << OGcartTree.prediction_errors(dt) << endl;

  cout << "Taille du cluster original " << cl.clusters.size() <<endl;
  int cpt = 0;
  while(keepGoing and clust.clusters.size() < dt.I){
    cpt += 1;
    dataset currentDt = clust.createDt(dt);

    Tree cartTree = pythonCART(currentDt, p.D,p.C,p.Nmin);
    
    solC.nbCl.push_back(clust.clusters.size());
    keepGoing = isIntersecting(cartTree, clust, dt, currentDt, true, true);

    cout << "Iteration " << cpt << endl;
    cout << "Taille du cluster " << clust.clusters.size() <<endl;
    cout << "Erreur sur Cl " << cartTree.prediction_errors(currentDt) << endl;
    cout << "Erreur sur Dt " <<cartTree.prediction_errors(dt) << endl;

    //cartTree.write_tree("iteration"+to_string(cpt)+".txt");

    solC.T.push_back(cartTree);
    solC.nbIter += 1;
    solC.errCl.push_back(solC.T[solC.T.size()-1].prediction_errors(currentDt));
    solC.errDt.push_back(solC.T[solC.T.size()-1].prediction_errors(dt));
  }

  solC.nbClFinal = clust.clusters.size();
  solC.finalCl = clust.clustering_copy();


  //clust.showClusters();

  return solC;
}

solClust iteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, double timeL){
  clustering clust = cl.clustering_copy();
  parameters pm = p.parameters_copy();

  pm.MIPFocus = 1;

  bool keepGoing = true;

  solClust solC = solClust();
  solC.nbClInitial = clust.clusters.size();

  //clust.showClusters();
  int cpt = 0;
  
  GRBEnv env = GRBEnv();
  while(timeL-solC.totTime> 10.0 and keepGoing and clust.clusters.size() < dt.I){
    cpt += 1;

    time_t t1 = time (NULL);
    
    dataset currentDt = clust.createDt(dt,true);
    pm.update(currentDt);
    pm.updateMu(dt); // pas pratique mais pas grave
    
    GRBModel md = GRBModel(env);
    build_model model = build_model(md,currentDt,clust,modelt,pm);


    // warmstart avec medoid en représentant c'est oké !!
    Tree cartTree = pythonCART(currentDt,pm.D,pm.C,pm.Nmin);
    if (cpt == 1){
      model.add_warmstart(cartTree, currentDt);
    }
    else{
      if (cartTree.prediction_errors(currentDt) > solC.T.back().prediction_errors(currentDt)){
	model.add_warmstart(solC.T.back(), currentDt);
      }
      else{
	model.add_warmstart(cartTree, currentDt);
      }
      }
    model.completeWarmStartForVarRep(clust);
    
    solution sol = model.solve2(md,timeL-solC.totTime);
    for (int i=0; i<dt.I; i++){
      if(model.var.r[i].get(GRB_DoubleAttr_X) > 0.99){
	clust.clusters[clust.placeOf[clust.clusterOf[i]]].rep = i;
      }
    }
    currentDt = clust.createDt(dt,false,true);
    
    time_t t2 = time (NULL);
    solC.nbCl.push_back(clust.clusters.size()); // on met le nombre de cl considérés avant de potentiellement changer les clusters
    keepGoing = isIntersecting(sol.T, clust, dt, currentDt, modelt.univ, pm.Nmin);
    time_t t3 = time (NULL);

    /*
    // on récupère des infos mais c'est à court terme
    freopen("yooooooooooooooooooooooo.txt", "a", stdout);
    cout << "Iteration " << cpt << endl;
    //cout << cartTree.prediction_errors(currentDt) << endl;
    cout << sol.obj << endl;
    cout << sol.gap << endl;
    cout << sol.T.prediction_errors(dt) << endl;
    cout << sol.time << endl;
    freopen("/dev/tty", "w", stdout);
    // a delete plus tard */

    solC.T.push_back(sol.T);
    solC.solveTime.push_back(t2-t1);
    solC.centeringTime.push_back(t3-t2);
    solC.gap.push_back(sol.gap);
    solC.totTime += t3 - t1;
    solC.nbIter += 1;
    solC.obj.push_back(sol.obj);
    solC.errCl.push_back(solC.T.back().prediction_errors(currentDt));
    solC.errDt.push_back(solC.T.back().prediction_errors(dt));

    if (solC.totTime > timeL){
      break;
    }
  }
  solC.nbClFinal = clust.clusters.size();
  solC.finalCl = clust.clustering_copy();

  if (solC.errCl.back() == 0){
    if (solC.errDt.back() == 0){
      solC.pseudoGap = 0.0;
    }
    else{
      solC.pseudoGap = 1.0;
    }
  }
  else{
    solC.pseudoGap = (double)(solC.errDt.back()-solC.errCl.back())/(double)solC.errDt.back();
  }
  
  return solC;
}

/*

solClust iteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, int timeL){ // A MODIFFFF
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
*/
