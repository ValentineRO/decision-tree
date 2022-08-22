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

training_results learning_Bertsimas(dataset& dt_train, dataset& dt_validation, dataset& dt_test, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin, double mu){

  training_results tr = training_results();

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
      parameters p = parameters(D,dt_train,0,mu,C,false,Nmin);

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
      
      tr.time += sol.time;
      tr.nb_opti += (int)(sol.time<= time_limit);
      tr.nb_iter += 1;
    }
  }

  int best = dom_trees.BestTree(true);
  tr.best_tree[0] = dom_trees.trees[best].namefile;
  tr.perc_err[0] = (double)dom_trees.trees[best].Et / (double) dt_test.I;
  tr.alph[0] = (dom_trees.alph[best]+dom_trees.alph[best+1])/2.0;
  tr.alph[0] /= (float)dt_train.L_h;

  return tr;
}

training_results learning(dataset& dt_train, dataset& dt_validation, dataset& dt_test, baseModel bm, bool univ, int DMAX, double time_limit, int Nmin, double mu){
  
  training_results tr = training_results();

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
      tr.nb_iter += 1;
      cout << "D = " << D << " et C = " << C << " / CMAX = " << CMAX << endl;
      freopen("gurobi_text.txt", "a", stdout);
      GRBEnv env = GRBEnv();
      parameters p = parameters(D,dt_train,1.0/(2*C),mu,C,false,Nmin);

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
      
      tr.time += sol.time;
      tr.nb_opti += (int)(sol.time<= time_limit);

      previousC = sol.nb_br;
      C = sol.nb_br-1;
    }
  }

  int* best = new int[4];
  
  best[0] = dom_trees_s_pp.BestTree(true);
  tr.best_tree[0] = dom_trees_s_pp.trees[best[0]].namefile;
  tr.perc_err[0] = (double) dom_trees_s_pp.trees[best[0]].Et / (double) dt_test.I;
  tr.alph[0] = (dom_trees_s_pp.alph[best[0]]+dom_trees_s_pp.alph[best[0]+1])/2.0;
  
  best[1] = dom_trees_s_pp.BestTree(false);
  tr.best_tree[1] = dom_trees_s_pp.trees[best[1]].namefile;
  tr.perc_err[1] = (double) dom_trees_s_pp.trees[best[1]].Et / (double) dt_test.I;
  tr.alph[1] = (dom_trees_s_pp.alph[best[1]]+dom_trees_s_pp.alph[best[1]+1])/2.0;
  
  best[2] = dom_trees_a_pp.BestTree(true);
  tr.best_tree[2] = dom_trees_a_pp.trees[best[2]].namefile;
  tr.perc_err[2] = (double) dom_trees_a_pp.trees[best[2]].Et / (double) dt_test.I;
  tr.alph[2] = (dom_trees_a_pp.alph[best[2]]+dom_trees_a_pp.alph[best[2]+1])/2.0;
  
  best[3] = dom_trees_a_pp.BestTree(false);
  tr.best_tree[3] = dom_trees_a_pp.trees[best[3]].namefile;
  tr.perc_err[3] = (double) dom_trees_a_pp.trees[best[3]].Et / (double) dt_test.I;
  tr.alph[3] = (dom_trees_a_pp.alph[best[3]]+dom_trees_a_pp.alph[best[3]+1])/2.0;
  
  // best between sans post-process et with post-process pour section du plus petit alpha
  if (dom_trees_s_pp.trees[best[0]].Ev < dom_trees_a_pp.trees[best[2]].Ev){
    tr.best_tree[4] = dom_trees_s_pp.trees[best[0]].namefile;
    tr.perc_err[4] = tr.perc_err[0];
    tr.alph[4] = tr.alph[0];
  }
  else{
    tr.best_tree[4] = dom_trees_a_pp.trees[best[2]].namefile;
    tr.perc_err[4] = tr.perc_err[2];
    tr.alph[4] = tr.alph[2];
  }
  // best between sans post-process et with post-process pour section du plus grand alpha
  if (dom_trees_s_pp.trees[best[1]].Ev < dom_trees_a_pp.trees[best[3]].Ev){
    tr.best_tree[5] = dom_trees_s_pp.trees[best[1]].namefile;
    tr.perc_err[5] = tr.perc_err[1];
    tr.alph[5] = tr.alph[1];
  }
  else{
    tr.best_tree[5] = dom_trees_a_pp.trees[best[3]].namefile;
    tr.perc_err[5] = tr.perc_err[3];
    tr.alph[5] = tr.alph[3];
  }

  for (int i=0; i<6; i++){
    tr.alph[i] /= (float)dt_train.L_h;
  }
  
  return tr;
}
