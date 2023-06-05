#include "train_tree.h"

string model_name(baseModel bm, bool univ){
  string s;
  if (bm == baseModel::F){s = "F";}
  if (bm == baseModel::OCT){s = "OCT";}
  if (bm == baseModel::QOCT){s = "QOCT";}
  if (bm == baseModel::FOCT){s = "FOCT";}
  if (bm == baseModel::QOCT){s = "GOCT";}
  if (!univ){s += "H";}
  return s;
}

training_results learning_Bertsimas(dataset& dt_train, dataset& dt_validation, dataset& dt_test, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin){

  training_results tr = training_results(DMAX-1);

  int nb_sol_tot = 0;
  
  for (int i = 1; i<= DMAX; i++){
    if (univ){
      nb_sol_tot += (int)pow(2,i)-i;
    }
    else{
      nb_sol_tot += (int)(pow(2,i)-1)*dt_train.J-i+1;
    }
  }

  model_type mt = model_type(bm, univ, false, true, false);

  dominating_trees dom_trees = dominating_trees(nb_sol_tot,dt_train.L_h);

  for (int D=1; D<=DMAX; D++){
    int CMAX = (int)pow(2,D)-1;
    if (!univ){ CMAX *= dt_train.J;}
    
    for (int C=D; C<=CMAX; C++){
      cout << "D = " << D << " et C = " << C << " / CMAX = " << CMAX << endl;
      freopen("gurobi_text.txt", "a", stdout);
      GRBEnv env = GRBEnv();
      parameters p = parameters(D,dt_train,0,C,false,Nmin);

      string CARTnamefile;
      if (univ){
	CARTnamefile = "../TreesAndPartitions/"+dt_train.name+"/CART_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(C)+"_Nmin.txt";
      }
      else{
	int CMAXuniv = min(C,CMAX/dt_train.J);
	CARTnamefile = "../TreesAndPartitions/"+dt_train.name+"/CART_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(CMAXuniv)+"_Nmin.txt";
      }

      Tree warmstartTree = Tree(CARTnamefile);
      int error_train_CART = warmstartTree.prediction_errors(dt_train);
      int best_tree = dom_trees.BestWarmstart(D,C);

      if (best_tree != -1){
	if (dom_trees.trees[best_tree].E < error_train_CART){
	  warmstartTree = Tree(dom_trees.trees[best_tree].namefile);
	}
      }
      
      GRBModel md = GRBModel(env);
      build_model model = build_model(md,dt_train,mt,p);
      model.add_warmstart(warmstartTree, dt_train);
      solution sol = model.solve(md,time_limit);

      string treeName = "../TreesAndPartitions/"+dt_train.name+"/sol_"+model_name(bm,univ)+"_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(C)+".txt";
      sol.T.write_tree(treeName);

      int error_val = sol.T.prediction_errors(dt_validation);
      int error_test = sol.T.prediction_errors(dt_test);
      
      optimal_tree opt = optimal_tree(D,sol.nb_br,(int)sol.obj,error_val,error_test,treeName);
      bool added = dom_trees.add_tree(opt);

      freopen("/dev/tty", "w", stdout);
      remove("gurobi_text.txt");
      
      tr.globalTime += sol.time;
      tr.globalNbOpti += (int)(sol.time<= time_limit);
      tr.globalNbIter += 1;
    }

    if (D != 1){
      tr.time[D-2] = tr.globalTime;
      tr.nbOpti[D-2] = tr.globalNbOpti;
      tr.nbIter[D-2] = tr.globalNbIter;
      
      int best = dom_trees.BestTree(true);      
      tr.bestTreeWithoutPP[D-2] = dom_trees.trees[best].namefile;
      tr.percErrWithoutPP[D-2] = (double)dom_trees.trees[best].Et / (double) dt_test.I;
      tr.alphWithoutPP[D-2] = (dom_trees.alph[best]+dom_trees.alph[best+1])/2.0;
      tr.alphWithoutPP[D-2] /= (float)dt_train.L_h;
    }
  }

  return tr;
}

training_results learning(dataset& dt_train, dataset& dt_validation, dataset& dt_test, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin){
  
  training_results tr = training_results(DMAX-1);

  int nb_sol_tot = 0;
  for (int i = 1; i<= DMAX; i++){
    if (univ){
      nb_sol_tot += (int)pow(2,i)-i;
    }
    else{
      nb_sol_tot += (int)pow(2,i)*dt_train.J-i;
    }
  }
  
  model_type mt = model_type(bm, univ, false, true, false);

  dominating_trees dom_trees_a_pp = dominating_trees(nb_sol_tot,dt_train.L_h),
    dom_trees_s_pp = dominating_trees(nb_sol_tot,dt_train.L_h);

  for (int D=1; D<=DMAX; D++){
    int CMAX = (int)pow(2,D)-1;
    if (!univ){ CMAX *= dt_train.J;}
    
    int C = CMAX;
    int previousC = -1;

    while (C >= D){
      tr.globalNbIter += 1;
      cout << "D = " << D << " et C = " << C << " / CMAX = " << CMAX << endl;
      freopen("gurobi_text.txt", "a", stdout);
      GRBEnv env = GRBEnv();
      parameters p = parameters(D,dt_train,1.0/(2*C),C,false,Nmin);

      string CARTnamefile;
      
      if (univ){
	CARTnamefile = "../TreesAndPartitions/"+dt_train.name+"/CART_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(C);
      }
      else{
	int CMAXuniv = min(C,CMAX/dt_train.J);
	CARTnamefile = "../TreesAndPartitions/"+dt_train.name+"/CART_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(CMAXuniv);
      }
      if (Nmin !=0){
	CARTnamefile = CARTnamefile + "_Nmin.txt";
      }
      else{
	CARTnamefile = CARTnamefile + ".txt";
      }
      Tree warmstartTree = Tree(CARTnamefile);
      int error_train_CART = warmstartTree.prediction_errors(dt_train);
      int best_tree = dom_trees_s_pp.BestWarmstart(D,C); // pas besoin de prendre la solution avec post-process (ça n'a pas d'importance)

      if (previousC != -1){ 
	Tree previousTree = Tree("../TreesAndPartitions/"+dt_train.name+"/sol_"+model_name(bm,univ)+"_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(previousC)+".txt");
	if (univ){
	  previousTree = previousTree.removeSplit(dt_train);
	}
	else{
	  previousTree = previousTree.reduceComplexity(dt_train,Nmin);
	}
	int errorModifiedTree = previousTree.prediction_errors(dt_train);
	if (errorModifiedTree < error_train_CART){
	  warmstartTree = previousTree;
	}
      }

      if (best_tree != -1){
	if (dom_trees_s_pp.trees[best_tree].E < error_train_CART){
	  warmstartTree = Tree(dom_trees_s_pp.trees[best_tree].namefile);
	}
      }

      GRBModel md = GRBModel(env);
      build_model model = build_model(md,dt_train,mt,p);
      model.add_warmstart(warmstartTree, dt_train);
      solution sol = model.solve(md,time_limit);

      string treeName = "../TreesAndPartitions/"+dt_train.name+"/sol_"+model_name(bm,univ)+"_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(sol.nb_br)+".txt";
      sol.T.write_tree(treeName);

      int error_train = (int)floor(sol.obj);

      // partie sans post-process
      int error_val = sol.T.prediction_errors(dt_validation);
      int error_test = sol.T.prediction_errors(dt_test);
      
      optimal_tree opt_s_pp = optimal_tree(D,sol.nb_br,error_train,error_val,error_test,treeName);
      bool added = dom_trees_s_pp.add_tree(opt_s_pp);

      // partie avec post-process

      if (univ){
	sol.T.post_processing_b(dt_train, bm!=baseModel::F); // les missclassif comptent si on a pas F en base model
      }
      else{
	sol.T.post_processing_a_b(dt_train, bm!=baseModel::F);
      }

      string postProcessedTreeName = "../TreesAndPartitions/"+dt_train.name+"/sol_"+model_name(bm,univ)+"_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(sol.nb_br)+"_PP.txt";
      sol.T.write_tree(postProcessedTreeName);

      error_val = sol.T.prediction_errors(dt_validation);
      error_test = sol.T.prediction_errors(dt_test);
      
      optimal_tree opt_a_pp = optimal_tree(D,sol.nb_br,error_train,error_val,error_test,postProcessedTreeName);
      added = dom_trees_a_pp.add_tree(opt_a_pp);

      freopen("/dev/tty", "w", stdout);
      remove("gurobi_text.txt");
      
      tr.globalTime += sol.time;
      tr.globalNbOpti += (int)(sol.time<= time_limit);

      previousC = sol.nb_br;
      C = sol.nb_br-1;
    }
    if (D != 1){
      tr.time[D-2] = tr.globalTime;
      tr.nbOpti[D-2] = tr.globalNbOpti;
      tr.nbIter[D-2] = tr.globalNbIter;

      int bestSansPP = dom_trees_s_pp.BestTree(true),
	bestAvecPP = dom_trees_a_pp.BestTree(true);

      tr.bestTreeWithoutPP[D-2] = dom_trees_s_pp.trees[bestSansPP].namefile;
      tr.percErrWithoutPP[D-2] = (double) dom_trees_s_pp.trees[bestSansPP].Et / (double) dt_test.I;
      tr.alphWithoutPP[D-2] = (dom_trees_s_pp.alph[bestSansPP]+dom_trees_s_pp.alph[bestSansPP+1])/2.0;
      tr.alphWithoutPP[D-2] /= (float)dt_train.L_h;

      tr.bestTreeWithPP[D-2] = dom_trees_a_pp.trees[bestAvecPP].namefile;
      tr.percErrWithPP[D-2] = (double) dom_trees_a_pp.trees[bestAvecPP].Et / (double) dt_test.I;
      tr.alphWithPP[D-2] = (dom_trees_a_pp.alph[bestAvecPP]+dom_trees_a_pp.alph[bestAvecPP+1])/2.0;
      tr.alphWithPP[D-2] /= (float)dt_train.L_h;

      if (dom_trees_s_pp.trees[bestSansPP].Ev < dom_trees_a_pp.trees[bestAvecPP].Ev){
	tr.bestTreeBestOfBoth[D-2] = tr.bestTreeWithoutPP[D-2];
	tr.percErrBestOfBoth[D-2] = tr.percErrWithoutPP[D-2];
	tr.alphBestOfBoth[D-2] = tr.alphWithoutPP[D-2];
      }
      else{
	tr.bestTreeBestOfBoth[D-2] = tr.bestTreeWithPP[D-2];
	tr.percErrBestOfBoth[D-2] = tr.percErrWithPP[D-2];
	tr.alphBestOfBoth[D-2] = tr.alphWithPP[D-2];
      }
    }
  }
  
  return tr;
}

training_results_cl learningWithClustering(dataset& dt_train, dataset& dt_validation, dataset& dt_test, clustering cl, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin, int selectingStrat){

  training_results_cl tr = training_results_cl(DMAX-1);

  int nb_sol_tot = 0;
  for (int i = 1; i<= DMAX; i++){
    if (univ){
      nb_sol_tot += (int)pow(2,i)-i;
    }
    else{
      nb_sol_tot += (int)pow(2,i)*dt_train.J-i;
    }
  }
  
  model_type mt = model_type(bm, univ, false, true, false);

  dominating_trees dom_trees = dominating_trees(nb_sol_tot,dt_train.L_h);
    
  for (int D=1; D<=DMAX; D++){
    int CMAX = (int)pow(2,D)-1;
    if (!univ){ CMAX *= dt_train.J;}
    
    int C = CMAX;
    int previousC = -1;

    while (C >= D){
      tr.globalNbIter += 1;
      cout << "D = " << D << " et C = " << C << " / CMAX = " << CMAX << endl;

      parameters p = parameters(D,dt_train,1.0/(2*C),C,false,Nmin);

      solClust solC = approxIteratingOTP(dt_train, cl, mt, p, time_limit);

      int bestTree = -1;
      if (selectingStrat == 0){// si pas opti (ie pseudoGap > 0), meilleur du train
	bestTree = solC.T.size()-1;
	if (solC.pseudoGap > 0){
	  for (int t=0; t<solC.T.size()-1; t++){
	    if (solC.errDt[bestTree] > solC.errDt[t]){
	      bestTree = t;
	    }
	  }
	}
      }
      if (selectingStrat == 1){// meilleur du train
	bestTree = solC.T.size()-1;
	for (int t=0; t<solC.T.size()-1; t++){
	  if (solC.errDt[bestTree] > solC.errDt[t]){
	    bestTree = t;
	  }
	}
      }
      if (selectingStrat == 2){// si pas opti (ie pseudoGap > 0), meilleur de la validation
	bestTree = solC.T.size()-1;
	if (solC.pseudoGap > 0){
	  int val = solC.T[bestTree].prediction_errors(dt_validation);
	  for (int t=0; t<solC.T.size()-1; t++){
	    int val_t = solC.T[t].prediction_errors(dt_validation);
	    if (val > val_t){
	      bestTree = t;
	      val = val_t;
	    }
	  }
	}

      }
      if (selectingStrat == 3){// meilleur de la validation
	bestTree = solC.T.size()-1;
	int val = solC.T[bestTree].prediction_errors(dt_validation);
	for (int t=0; t<solC.T.size()-1; t++){
	  int val_t = solC.T[t].prediction_errors(dt_validation);
	  if (val > val_t){
	    bestTree = t;
	    val = val_t;
	  }
	}
      }

      int nb_br = solC.T[bestTree].countComplexity();
      
      string itOptName = "../TreesAndPartitions/"+dt_train.name+"/itOpt_"+model_name(bm,univ)+"_"+cl.name+"_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(nb_br)+"_selStrat"+to_string(selectingStrat)+".txt";
      solC.write(itOptName);

      string treeName = "../TreesAndPartitions/"+dt_train.name+"/sol_"+model_name(bm,univ)+"_"+cl.name+"_part"+to_string(dt_train.partition)+"_D"+to_string(D)+"_C"+to_string(nb_br)+"_selStrat"+to_string(selectingStrat)+".txt";
      solC.T[bestTree].write_tree(treeName);

      int error_train = (int)floor(solC.errDt[bestTree]);
      int error_val = solC.T[bestTree].prediction_errors(dt_validation);
      int error_test = solC.T[bestTree].prediction_errors(dt_test);

      optimal_tree opt = optimal_tree(D,nb_br,error_train,error_val,error_test,treeName);
      bool added = dom_trees.add_tree(opt);
      
      tr.globalTime += solC.totTime;
      tr.globalNbOpti += (int)(solC.totTime<= time_limit);
      tr.globalNbIterCl += solC.nbIter;
      tr.globalPseudoGap += solC.pseudoGap;
      tr.globalFinalReduction += (double)solC.nbClFinal/(double)dt_train.I;

      previousC = nb_br;
      C = nb_br-1;
    }
    if (D != 1){
      tr.time[D-2] = tr.globalTime;
      tr.nbOpti[D-2] = tr.globalNbOpti;
      tr.nbIter[D-2] = tr.globalNbIter;
      tr.nbIterCl[D-2] = (double)tr.globalNbIterCl / (double)tr.globalNbIter;
      tr.pseudoGap[D-2] = (double)tr.globalPseudoGap / (double)tr.globalNbIter;
      tr.finalReduction[D-2] = (double)tr.globalFinalReduction / (double)tr.globalNbIter;

      int bestTree = dom_trees.BestTree(true);

      tr.bestTree[D-2] = dom_trees.trees[bestTree].namefile;
      tr.errTrain[D-2] = dom_trees.trees[bestTree].E;
      tr.errTest[D-2] = dom_trees.trees[bestTree].Et;
    }
  }
  
  return tr;
}
