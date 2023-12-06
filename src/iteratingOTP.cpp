#include "iteratingOTP.h"

extern GurobiEnvironment& gurobiEnv;

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

bool isIntersectingBasic(Tree T, clustering& cl, dataset& initialDt, dataset& currentDt, bool variableRep){
  // currentDt must be made up of representatives

  int* predData = new int[initialDt.I];
  T.predict_leaves(initialDt,predData);
  //int* predRep = new int[currentDt.I];

  bool isIntersecting = false;
  for (int g=0; g<currentDt.I; g++){
    vector<vector<int>> div;
    vector<int> val;

    //div.push_back({});
    //val.push_back(predRep[g]);

    for (auto i : cl.clusters[g].pts){
      auto it = find(val.begin(), val.end(), predData[i]);
      if (it != val.end()){
	int index = distance(val.begin(), it);
	div[index].push_back(i);
      }
      else{
	div.push_back({i});
	val.push_back(predData[i]);
      }
    }

    if (div.size() > 1){
      cl.breakCluster(g, div, initialDt, variableRep);
      // the following lines are added because it is not done in break clusters because of other isIntersecting functions
      cl.clusters[g].computeBarycenter(initialDt);
      if (variableRep){
	cl.clusters[g].rep = cl.clusters[g].m;
      }
      isIntersecting = true;
    }
  }
  return isIntersecting;
}

struct RealNumberCompare { // used in isIntersectingUNIV
  bool operator() (const double& a, const double& b) const {
    return a < b;
  }
};

double roundToPrecision(double value, double precision) { // used in isIntersectingUNIV
    return round(value / precision) * precision;
}

bool isIntersectingUNIV(Tree T, clustering& cl, dataset& initialDt, dataset& currentDt, bool variableRep){
  bool* goodPred = new bool[currentDt.I];
  T.isGoodPrediction(currentDt, goodPred);

  /*
  for (int g=0; g<currentDt.I; g++){
    cout << "Cluster " << g << " - " << goodPred[g] <<endl;
    }*/
  
  bool B = false;

  for (int t=0; t<T.N; t++){
    vector<int> IL = {},
      IR = {};
    T.compute_ILIR_t(IL, IR, t, currentDt, true);
    
    if (T.c[t] == -1 and IL.size()!=0 and IR.size()!=0){
      vector<int> leftPoints = {},
	rightPoints = {};

      map<int,bool> canBeCut;
      
      // we identify points we want to send left and right
      for (auto g: IL){
	if (goodPred[g]){
	  canBeCut[g] = true;
	  for (auto i: cl.clusters[g].pts){
	    leftPoints.push_back(i);
	  }
	}
	else{
	  canBeCut[g] = false;
	  for (auto i: cl.clusters[g].pts){
	    bool goodPredIfSw = T.goodPredIfSwap(initialDt, i, t);
	    if (goodPredIfSw){
	      canBeCut[g] = true;
	      rightPoints.push_back(i);
	    }
	  }
	}
	if (cl.clusters[g].pts.size() == 1){
	  canBeCut[g] = false;
	}
      }
      for (auto g: IR){
	if (goodPred[g]){
	  canBeCut[g] = true;
	  for (auto i: cl.clusters[g].pts){
	    rightPoints.push_back(i);
	  }
	}
	else{
	  canBeCut[g] = false;
	  for (auto i: cl.clusters[g].pts){
	    bool goodPredIfSw = T.goodPredIfSwap(initialDt, i, t);
	    if (goodPredIfSw){
	      canBeCut[g] = true;
	      leftPoints.push_back(i);
	    }
	  }
	}
	if (cl.clusters[g].pts.size() == 1){
	  canBeCut[g] = false;
	}
      }

      int sepAtt = -1; // identify the attribute on which we separate data points
      for (int j=0; j<T.J; j++){
	if (T.a[t*T.J+j] == 1.0){
	  sepAtt = j;
	  break;
	}
      }
      double maxLeft = 0;
      for (auto i: leftPoints){
	maxLeft = max(maxLeft, roundToPrecision(initialDt.X[i*T.J+sepAtt],mu/10));
      }


      double minRight = 1;
      for (auto i: rightPoints){
	minRight = min(minRight, roundToPrecision(initialDt.X[i*T.J+sepAtt],mu/10));
      }

      //cout << "Noeud " << t << " - minRight " << minRight << " - maxLeft " << maxLeft<< endl;

      if (minRight > maxLeft){
	T.b[t] = (minRight+maxLeft)/2;
      }
      else{
	// we want to know where representatives go 
	double maxRepLeft = 0;
	for (auto g: IL){
	  if (goodPred[g]){
	    if (variableRep){
	      maxRepLeft = max(maxRepLeft, roundToPrecision(initialDt.X[cl.clusters[g].rep*T.J+sepAtt], mu/10));
	    }
	    else{
	      maxRepLeft = max(maxRepLeft, roundToPrecision(currentDt.X[g*T.J+sepAtt], mu/10));
	    }
	  }
	}
	double minRepRight = 1;
	for (auto g: IR){
	  if (goodPred[g]){
	    if (variableRep){
	      minRepRight = min(minRepRight, roundToPrecision(initialDt.X[cl.clusters[g].rep*T.J+sepAtt], mu/10));
	    }
	    else{
	      minRepRight = min(minRepRight, roundToPrecision(currentDt.X[g*T.J+sepAtt], mu/10)); 
	    }
	  }
	}
	
	// we can know check which values of X are critical
	map<double,pair<int,int>, RealNumberCompare> Xval;

	int cantFollow = 0;
	
	for (auto i: leftPoints){
	  double val = roundToPrecision(initialDt.X[i*T.J+sepAtt],mu/10);
	  if (val > maxRepLeft){
	    if (val < minRepRight){
	      Xval[val].first += 1;
	    }
	    else{
	      cantFollow += 1;
	    }
	  }
	}
	
	int notFollowing = 0;  // we start computing the number of data not following with the left most value possible
	for (auto i: rightPoints){
	  double val = roundToPrecision(initialDt.X[i*T.J+sepAtt],mu/10);
	  if (val < minRepRight){
	    if (val > maxRepLeft){
	      Xval[val].second += 1;
	      notFollowing += 1;
	    }
	    else{
	      cantFollow += 1;
	    }
	  }
	}

	double prevVal = maxRepLeft;
	
	double bestInterv = 0;
	double bestNotFoll = initialDt.I;
	double sepVal;
	double interv;

	for (const auto& it : Xval) {
	  if (notFollowing < bestNotFoll){
	    interv = it.first - prevVal;
	    if (bestInterv < interv){
	      bestInterv = interv;
	      bestNotFoll = notFollowing;
	      sepVal = (it.first + prevVal)/2;
	    }
	  }
	  // update for next iteration
	  prevVal = it.first;
	  notFollowing += it.second.first; // now that we go to next biggest value, elements that were suppose to go left are not following their representatives
	  notFollowing -= it.second.second; // now that we go to next biggest value, elements that were suppose to go right are following their representatives
	}
	if (notFollowing < bestNotFoll){ // this is the last interval
	  interv = minRepRight - prevVal;
	  if (bestInterv < interv){
	    bestInterv = interv;
	    bestNotFoll = notFollowing;
	    sepVal = (minRepRight + prevVal)/2;
	  }
	}

	for (auto g: IL){
	  if (canBeCut[g]){
	    bool refGoingLeft = true;
	    if(!goodPred[g]){
	      double sum_ax = 0;
	      for (int j=0; j<initialDt.J; j++){
		sum_ax += currentDt.X[g*currentDt.J+j]*T.a[t*currentDt.J+j];
	      }

	      refGoingLeft = lessThan(sum_ax, T.b[t]);
	    }
	    vector<vector<int>> division;
	    division.push_back({});
	    division.push_back({});
	    for (auto i: cl.clusters[g].pts){
	      double sum_ax = 0;
	      for (int j=0; j<initialDt.J; j++){
		sum_ax += initialDt.X[i*initialDt.J+j]*T.a[t*initialDt.J+j];
	      }
	      if (lessThan(sum_ax, T.b[t])){
		if (refGoingLeft){
		  division[0].push_back(i);
		}
		else{
		  division[1].push_back(i);
		}
	      }
	      else{
		if (refGoingLeft){
		  division[1].push_back(i);
		}
		else{
		  division[0].push_back(i);
		}
	      }
	    }

	    if (division[0].size()>0 and division[1].size()>0){
	      cl.breakCluster(g, division, initialDt, variableRep);
	      B = true; // there is a cluster break !!
	    }
	    else{
	      if (division[0].size()==0){  // can only happend if !variableRep (and some data points have already left the clusters)
		cl.clusters[g].computeBarycenter(initialDt); // in this case we modify the representative
		cl.clusters[g].computeMedoid(initialDt);
		for (int j=0; j<initialDt.J; j++){
		  currentDt.X[g*currentDt.J+j] = cl.clusters[g].bar[j];
		}
		if (!B){
		  cout << "wahou c pas du tout normal" << endl; // just in case we print an error message

		}
	      }
	    }
	    
	  }
	  
	}
	for (auto g: IR){
	  if (canBeCut[g]){
	    bool refGoingRight = true;
	    if(!goodPred[g]){
	      double sum_ax = 0;
	      for (int j=0; j<currentDt.J; j++){
		sum_ax += currentDt.X[g*currentDt.J+j]*T.a[t*currentDt.J+j];
	      }

	      refGoingRight = greaterThanOrEqualTo(sum_ax, T.b[t]);
	    }
	    vector<vector<int>> division;
	    division.push_back({});
	    division.push_back({});
	    for (auto i: cl.clusters[g].pts){
	      double sum_ax = 0;
	      for (int j=0; j<initialDt.J; j++){
		sum_ax += initialDt.X[i*initialDt.J+j]*T.a[t*initialDt.J+j];
	      }
	      if (greaterThanOrEqualTo(sum_ax, T.b[t])){
		if (refGoingRight){
		  division[0].push_back(i);
		}
		else{
		  division[1].push_back(i);
		}
	      }
	      else{
		if (refGoingRight){
		  division[1].push_back(i);
		}
		else{
		  division[0].push_back(i);
		}
	      }
	    }
	    if (division[0].size()>0 and division[1].size()>0){
	      cl.breakCluster(g, division, initialDt, variableRep);
	      B = true; // there is a cluster break !!
	    }
	    else{
	      if (division[0].size()==0){  // can only happend if !variableRep (and some data points have already left the clusters)
		cl.clusters[g].computeBarycenter(initialDt); // in this case we modify the representative
		cl.clusters[g].computeMedoid(initialDt);
		for (int j=0; j<initialDt.J; j++){
		  currentDt.X[g*currentDt.J+j] = cl.clusters[g].bar[j];
		}
		if (!B){
		  cout << "wahou c pas du tout normal" << endl; // just in case we print an error message
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  for (int c=0; c<currentDt.I; c++){ // the original clusters
    cl.clusters[c].computeBarycenter(initialDt);
    cl.clusters[c].computeMedoid(initialDt);
  }

  return B;
}

bool isIntersecting(Tree T, clustering& cl, dataset& initialDt, dataset& currentDt, bool univ, int Nmin, bool variableRep){
  bool* goodPred = new bool[currentDt.I];
  T.isGoodPrediction(currentDt, goodPred);
  
  bool B = false;
  //GRBEnv env = GRBEnv();
  GRBEnv& env = gurobiEnv.getEnvironment();
  for (int t=0; t<T.N; t++){
    vector<int> IL = {},
      IR = {};
    T.compute_ILIR_t(IL, IR, t, currentDt, true);

    //cout << "noeud " << t << " - donnée à g: " << IL.size() << " - donnée à d: " << IR.size() << endl;
    
    if (T.c[t] == -1 and IL.size()!=0 and IR.size()!=0){
      vector<int> leftPoints = {},
	rightPoints = {};

      map<int,bool> canBeCut;
      
      // we identify points we want to send left and right
      int nbPoints = 0;
      for (auto g: IL){
	if (goodPred[g]){
	  canBeCut[g] = true;
	  for (auto i: cl.clusters[g].pts){
	    leftPoints.push_back(i);
	  }
	  nbPoints += cl.clusters[g].pts.size();
	}
	else{
	  canBeCut[g] = false;
	  for (auto i: cl.clusters[g].pts){
	    bool goodPredIfSw = T.goodPredIfSwap(initialDt, i, t);
	    if (goodPredIfSw){
	      canBeCut[g] = true;
	      rightPoints.push_back(i);
	      nbPoints += 1;
	    }
	  }
	}
	if (cl.clusters[g].pts.size() == 1){
	  canBeCut[g] = false;
	}
      }
      for (auto g: IR){
	if (goodPred[g]){
	  canBeCut[g] = true;
	  for (auto i: cl.clusters[g].pts){
	    rightPoints.push_back(i);
	  }
	  nbPoints += cl.clusters[g].pts.size();
	}
	else{
	  canBeCut[g] = false;
	  for (auto i: cl.clusters[g].pts){
	    bool goodPredIfSw = T.goodPredIfSwap(initialDt, i, t);
	    if (goodPredIfSw){
	      canBeCut[g] = true;
	      leftPoints.push_back(i);
	      nbPoints += 1;
	    }
	  }
	}
	if (cl.clusters[g].pts.size() == 1){
	  canBeCut[g] = false;
	}
      }

      GRBModel m = GRBModel(env);

      GRBVar* z_var;
      z_var = m.addVars(nbPoints, GRB_BINARY);

      GRBVar b1_var, b2_var;
      b1_var = m.addVar(-1, 1.0, 0.0, GRB_CONTINUOUS, "b1");
      b2_var = m.addVar(-1, 1.0, 0.0, GRB_CONTINUOUS, "b2");
      m.addConstr(b1_var,GRB_LESS_EQUAL,b2_var-mu,"b1<=b2-mu");
	
      GRBVar* a_var;
      GRBVar* ah_var;
      
      if (univ){
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

	for (int j=0; j<initialDt.J; j++){
	  m.addConstr(a_var[j],GRB_EQUAL,T.a[t*initialDt.J+j],"a_j=a*_j_j="+to_string(j));
	}
      }
      else{
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

	 int* sj_star = new int[initialDt.J];
	 for (int j=0; j<initialDt.J; j++){
	   if (equalTo(T.a[t*initialDt.J+j],0)){
	     sj_star[j] = 0;
	   }
	   else{
	     sj_star[j] = 1;
	   }
	 }

	 GRBLinExpr sum_ahj;
	 for (int j=0; j<T.J; j++){
	   m.addConstr(-sj_star[j],GRB_LESS_EQUAL,a_var[j],"-sj<=aj_j="+to_string(j));
	   m.addConstr(sj_star[j],GRB_GREATER_EQUAL,a_var[j],"sj>=aj_j="+to_string(j));
	   m.addConstr(-ah_var[j],GRB_LESS_EQUAL,a_var[j],"-ah_j<=aj_j="+to_string(j));
	   m.addConstr(ah_var[j],GRB_GREATER_EQUAL,a_var[j],"-ah_j>=aj_j="+to_string(j));
	   
	   sum_ahj += ah_var[j];
	 }
	 m.addConstr(sum_ahj,GRB_LESS_EQUAL, 1, "sum_ahj<=1");
      }

      string constraint_name;
      if (!variableRep){ // path of well predicted representatives (in case where rep is not a point of the cluster)
	for (auto g: IL){
	  if (goodPred[g]){
	    GRBLinExpr sum_ax = 0;
	    for (int j=0; j<initialDt.J; j++){
	      sum_ax += a_var[j]*cl.clusters[g].bar[j];
	    }
        
	    constraint_name = "representative_g_goes_into_left_branch_g="+to_string(g);
	    m.addConstr(sum_ax,GRB_LESS_EQUAL,b1_var,constraint_name);
	  }
	}

	for (auto g: IR){
	  if (goodPred[g]){
	    GRBLinExpr sum_ax = 0;
	    for (int j=0; j<initialDt.J; j++){
	      sum_ax += a_var[j]*cl.clusters[g].bar[j];
	    }
        
	    constraint_name = "representative_g_goes_into_right_branch_g="+to_string(g);
	    m.addConstr(sum_ax,GRB_GREATER_EQUAL,b2_var,constraint_name);
	  }
	}
      }
      else{
	for (auto g: IL){
	  if (goodPred[g]){
	    GRBLinExpr sum_ax = 0;
	    for (int j=0; j<initialDt.J; j++){
	      sum_ax += a_var[j]*initialDt.X[cl.clusters[g].rep*initialDt.J+j];
	    }
        
	    constraint_name = "representative_g_goes_into_left_branch_g="+to_string(g);
	    m.addConstr(sum_ax,GRB_LESS_EQUAL,b1_var,constraint_name);
	  }
	}

	for (auto g: IR){
	  if (goodPred[g]){
	    GRBLinExpr sum_ax = 0;
	    for (int j=0; j<initialDt.J; j++){
	      sum_ax += a_var[j]*initialDt.X[cl.clusters[g].rep*initialDt.J+j];
	    }
        
	    constraint_name = "representative_g_goes_into_right_branch_g="+to_string(g);
	    m.addConstr(sum_ax,GRB_GREATER_EQUAL,b2_var,constraint_name);
	  }
	}
      }

      double bigMLeft = 1+mu,
	bigMRight = 1;
      if (!univ){
	bigMLeft += 1;
	bigMRight += 1;
      }
      
      GRBLinExpr objFunc = 0;
      
      int cpt = 0;
      for (auto i: leftPoints){
	GRBLinExpr sum_ax = 0;
	for (int j=0; j<initialDt.J; j++){
	  sum_ax += a_var[j]*initialDt.X[i*initialDt.J+j];
	}
	
	constraint_name = "point_i_goes_into_left_branch_i="+to_string(i);
	m.addConstr(sum_ax, GRB_LESS_EQUAL, b1_var + bigMLeft*(1-z_var[cpt]),constraint_name);
	constraint_name = "point_i_goes_into_right_branch_i="+to_string(i);
	m.addConstr(sum_ax, GRB_GREATER_EQUAL, b2_var - bigMRight*z_var[cpt],constraint_name);

	objFunc += 1 - z_var[cpt];
	
	cpt++;
      }
      for (auto i: rightPoints){
	GRBLinExpr sum_ax = 0;
	for (int j=0; j<initialDt.J; j++){
	  sum_ax += a_var[j]*initialDt.X[i*initialDt.J+j];
	}

	constraint_name = "point_i_goes_into_left_branch_i="+to_string(i);
	m.addConstr(sum_ax, GRB_LESS_EQUAL, b1_var + bigMLeft*(1-z_var[cpt]),constraint_name);
	constraint_name = "point_i_goes_into_right_branch_i="+to_string(i);
	m.addConstr(sum_ax, GRB_GREATER_EQUAL, b2_var - bigMRight*z_var[cpt],constraint_name);

	objFunc += z_var[cpt];
	
	cpt++;
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

      T.b[t] = (b1_var.get(GRB_DoubleAttr_X) + b2_var.get(GRB_DoubleAttr_X))/ 2;

      if (univ){
	for (int j=0; j<initialDt.J; j++){
	  T.a[t*initialDt.J+j] = abs(a_var[j].get(GRB_DoubleAttr_X));
	}
      }
      else{
	for (int j=0; j<initialDt.J; j++){
	  T.a[t*initialDt.J+j] = a_var[j].get(GRB_DoubleAttr_X);
	}
      }
	
      if (greaterThan(m.get(GRB_DoubleAttr_ObjVal),0)){
	for (auto g: IL){
	  if (canBeCut[g]){
	    bool refGoingLeft = true;
	    if(!goodPred[g]){
	      double sum_ax = 0;
	      for (int j=0; j<initialDt.J; j++){
		sum_ax += currentDt.X[g*currentDt.J+j]*T.a[t*currentDt.J+j];
	      }

	      refGoingLeft = lessThan(sum_ax, T.b[t]);
	    }
	    vector<vector<int>> division;
	    division.push_back({});
	    division.push_back({});
	    for (auto i: cl.clusters[g].pts){
	      double sum_ax = 0;
	      for (int j=0; j<initialDt.J; j++){
		sum_ax += initialDt.X[i*initialDt.J+j]*T.a[t*initialDt.J+j];
	      }
	      if (lessThan(sum_ax, T.b[t])){
		if (refGoingLeft){
		  division[0].push_back(i);
		}
		else{
		  division[1].push_back(i);
		}
	      }
	      else{
		if (refGoingLeft){
		  division[1].push_back(i);
		}
		else{
		  division[0].push_back(i);
		}
	      }
	    }
	    if (division[0].size()>0 and division[1].size()>0){
	      cl.breakCluster(g, division, initialDt, variableRep);
	      B = true; // there is a cluster break !!
	    }
	    else{
	      if (division[0].size()==0){  // can only happend if !variableRep (and some data points have already left the clusters)
		cl.clusters[g].computeBarycenter(initialDt); // in this case we modify the representative
		cl.clusters[g].computeMedoid(initialDt);
		for (int j=0; j<initialDt.J; j++){
		  currentDt.X[g*currentDt.J+j] = cl.clusters[g].bar[j];
		}
		if (!B){
		  cout << "wahou c pas du tout normal" << endl; // just in case we print an error message
		}
	      }
	    }
	    
	  }
	}
	for (auto g: IR){
	  if (canBeCut[g]){
	    bool refGoingRight = true;
	    if(!goodPred[g]){
	      double sum_ax = 0;
	      for (int j=0; j<currentDt.J; j++){
		sum_ax += currentDt.X[g*currentDt.J+j]*T.a[t*currentDt.J+j];
	      }

	      refGoingRight = greaterThanOrEqualTo(sum_ax, T.b[t]);
	    }
	    vector<vector<int>> division;
	    division.push_back({});
	    division.push_back({});
	    for (auto i: cl.clusters[g].pts){
	      double sum_ax = 0;
	      for (int j=0; j<initialDt.J; j++){
		sum_ax += initialDt.X[i*initialDt.J+j]*T.a[t*initialDt.J+j];
	      }
	      if (greaterThanOrEqualTo(sum_ax, T.b[t])){
		if (refGoingRight){
		  division[0].push_back(i);
		}
		else{
		  division[1].push_back(i);
		}
	      }
	      else{
		if (refGoingRight){
		  division[1].push_back(i);
		}
		else{
		  division[0].push_back(i);
		}
	      }
	    }
	    if (division[0].size()>0 and division[1].size()>0){
	      cl.breakCluster(g, division, initialDt, variableRep);
	      B = true; // there is a cluster break !!
	    }
	    else{
	      if (division[0].size()==0){  // can only happend if !variableRep (and some data points have already left the clusters)
		cl.clusters[g].computeBarycenter(initialDt); // in this case we modify the representative
		cl.clusters[g].computeMedoid(initialDt);
		for (int j=0; j<initialDt.J; j++){
		  currentDt.X[g*currentDt.J+j] = cl.clusters[g].bar[j];
		}
		if (!B){
		  cout << "wahou c pas du tout normal" << endl; // just in case we print an error message
		}
	      }
	    }
	  }
	}
      }
    }
  }

  for (int c=0; c<currentDt.I; c++){ // the original clusters
    cl.clusters[c].computeBarycenter(initialDt);
    cl.clusters[c].computeMedoid(initialDt);
  }
  
  return B;
}


solClust approxIteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, bool noSplitChange, double timeL){
  clustering clust = cl.clustering_copy();
  parameters pm = p.parameters_copy();

  pm.MIPFocus = 1;

  bool keepGoing = true;

  solClust solC = solClust();
  solC.nbClInitial = clust.clusters.size();

  //clust.showClusters();
  int cpt = 0;
  
  //GRBEnv env = GRBEnv();
  GRBEnv& env = gurobiEnv.getEnvironment();
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
    /*
    clust.showClusters();
    dt.writeDataset("initialD.txt");
    currentDt.writeDataset("currD.txt");*/
    
    solution sol = model.solve2(md,timeL-solC.totTime);
    
    time_t t2 = time (NULL);
    solC.nbCl.push_back(clust.clusters.size()); // on met le nombre de cl considérés avant de potentiellement changer les clusters
    if (noSplitChange){
      keepGoing = isIntersectingBasic(sol.T, clust, dt, currentDt, false);
    }
    else{
      if (modelt.univ){
	keepGoing = isIntersectingUNIV(sol.T, clust, dt, currentDt);
      }
      else{
	keepGoing = isIntersecting(sol.T, clust, dt, currentDt, modelt.univ, pm.Nmin);
      }
      //keepGoing = isIntersecting(sol.T, clust, dt, currentDt, modelt.univ, pm.Nmin);
    }
    time_t t3 = time (NULL);
    
    time_t prevTime12 = t2-t1,
      prevTime23 = t3-t2;

    if (sol.gap > 0 and !keepGoing){ // we want to keep advancing if there is no cut but the optimality was not prove for the iteration
      solC.totTime += t3 - t1;
      time_t t1 = time (NULL);

      sol = model.solve2(md,timeL-solC.totTime,1.0);

      //sol.T.showTree();

      time_t t2 = time (NULL);
      solC.nbCl.back() = clust.clusters.size();

      if (noSplitChange){
      keepGoing = isIntersectingBasic(sol.T, clust, dt, currentDt, false);
      }
      else{
	if (modelt.univ){
	  keepGoing = isIntersectingUNIV(sol.T, clust, dt, currentDt);
	}
	else{
	  keepGoing = isIntersecting(sol.T, clust, dt, currentDt, modelt.univ, pm.Nmin);
	}
      }
    }
    else{
      prevTime12 = 0;
      prevTime23 = 0;
    }

    /*
    // on récupère des infos mais c'est à court terme
    //freopen("yoooooooooooo.txt", "a", stdout);
    cout << "Iteration " << cpt << endl;
    cout << cartTree.prediction_errors(currentDt) << endl;
    cout << sol.obj << endl;
    cout << sol.gap << endl;
    cout << sol.T.prediction_errors(dt) << endl;
    cout << sol.time << endl;
    //freopen("/dev/tty", "w", stdout);
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

    if (solC.totTime > timeL or solC.errDt.back() == 0){
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
    keepGoing = isIntersectingUNIV(cartTree, clust, dt, currentDt, false);

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

solClust iteratingOTP(dataset& dt, clustering& cl, model_type modelt, parameters p, bool noSplitChange, double timeL){
  clustering clust = cl.clustering_copy();
  parameters pm = p.parameters_copy();

  pm.MIPFocus = 1;

  bool keepGoing = true;

  solClust solC = solClust();
  solC.nbClInitial = clust.clusters.size();

  //clust.showClusters();
  int cpt = 0;
  
  //GRBEnv env = GRBEnv();
  GRBEnv& env = gurobiEnv.getEnvironment();
  while(timeL-solC.totTime> 10.0 and keepGoing and clust.clusters.size() < dt.I){
    cpt += 1;
    //cout << "Iteration " << cpt << endl;

    time_t t1 = time (NULL);
    
    dataset currentDt = clust.createDt(dt,true);
    pm.update(currentDt);
    pm.updateMu(dt); // pas pratique mais pas grave
    
    GRBModel md = GRBModel(env);
    build_model model = build_model(md,currentDt,clust,modelt,pm);

    //clust.showClusters();

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

    //sol.T.showTree();

    time_t t2 = time (NULL);
    solC.nbCl.push_back(clust.clusters.size()); // on met le nombre de cl considérés avant de potentiellement changer les clusters

    if (noSplitChange){
      keepGoing = isIntersectingBasic(sol.T, clust, dt, currentDt, true);
    }
    else{
      if (modelt.univ){
	keepGoing = isIntersectingUNIV(sol.T, clust, dt, currentDt, true);
      }
      else{
	keepGoing = isIntersecting(sol.T, clust, dt, currentDt, modelt.univ, pm.Nmin, true);
      }
    }

    
    time_t t3 = time (NULL);

    time_t prevTime12 = t2-t1,
      prevTime23 = t3-t2;

    if (sol.gap > 0 and !keepGoing and timeL > 10+solC.totTime+t3-t1){ // we want to keep advancing if there is no cut but the optimality was not prove for the iteration
      solC.totTime += t3 - t1;
      time_t t1 = time (NULL);

      sol = model.solve2(md,timeL-solC.totTime,1.0);
      for (int i=0; i<dt.I; i++){
	if(model.var.r[i].get(GRB_DoubleAttr_X) > 0.99){
	  clust.clusters[clust.placeOf[clust.clusterOf[i]]].rep = i;
	}
      }

      currentDt = clust.createDt(dt,false,true);

      //sol.T.showTree();

      time_t t2 = time (NULL);
      solC.nbCl.back() = clust.clusters.size();

      if (noSplitChange){
	keepGoing = isIntersectingBasic(sol.T, clust, dt, currentDt, true);
      }
      else{
	if (modelt.univ){
	  keepGoing = isIntersectingUNIV(sol.T, clust, dt, currentDt, true);
	}
	else{
	  keepGoing = isIntersecting(sol.T, clust, dt, currentDt, modelt.univ, pm.Nmin, true);
	}
      }

      time_t t3 = time (NULL);
    }
    else{
      prevTime12 = 0;
      prevTime23 = 0;
    }

    //sol.T.showTree();
    // on récupère des infos mais c'est à court terme

    /*
    //freopen("yooooooooooooooooooooooo.txt", "a", stdout);
    cout << "Iteration " << cpt << endl;
    cout << clust.clusters.size() << endl;
    cout << sol.T.prediction_errors(currentDt) << endl;
    cout << sol.obj << endl;
    cout << sol.gap << endl;
    cout << sol.T.prediction_errors(dt) << endl;
    cout << sol.time << endl;
    cout << t3 - t2 << endl;
    //freopen("/dev/tty", "w", stdout);
    */
    // a delete plus tard 

    solC.T.push_back(sol.T);
    solC.solveTime.push_back(prevTime12 + t2 - t1);
    solC.centeringTime.push_back(prevTime23 + t3 - t2);
    solC.gap.push_back(sol.gap);
    solC.totTime += prevTime12 + prevTime23 + t3 - t1;
    solC.nbIter += 1;
    solC.obj.push_back(sol.obj);
    solC.errCl.push_back((int)sol.obj);
    solC.errDt.push_back(solC.T.back().prediction_errors(dt));

    if (solC.totTime > timeL or solC.errDt.back()==0){
      break;
    }

    if (sol.gap < 0.0001 and solC.errCl.back()==solC.errDt.back()){
      break;
    }
  }

  //clust.showClusters();
  
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


