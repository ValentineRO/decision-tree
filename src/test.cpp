#include "test.h"

info_matrix::info_matrix(int n, int l, string col_names[]){
  nb_columns = n;
  nb_lines = l;
  columns_name = col_names;
  content = new string[nb_columns*nb_lines];
}

void info_matrix::write_line(int line_num, string line[]){
  for (int c=0; c<nb_columns; c++){
    content[line_num*nb_columns+c] = line[c];
  }
}

void info_matrix::write_csv(string namefile, int beg_line, int end_l){
  int end_line = nb_lines;
  if (end_l !=-1){
    end_line = end_l;
  }
  
  fstream file;
  
  file.open(namefile, ios::out);
  
  for (int c = 0; c < nb_columns; c++) {
    file << columns_name[c] << ",";
  }
  file << "\n";
  
  for (int l=beg_line ; l < end_line; l++) {
    for (int c = 0; c < nb_columns; c++) {
      file << content[l * nb_columns + c] << ",";
    }
    file << "\n";
  }
  file.close();
}


statistical_test_matrix::statistical_test_matrix(int nb_mod, int nb_p, int nb_dep, string dts, int Nm, float alph, float tl){
  nb_columns = 16;
  columns_name = new string[16];
  string cl_n[16] = { "Dataset", "Partition", "Model", "Time_Limit", "Profondeur", "Alpha", "Nmin", "Objective", "Error train", "Number of branchings", "Time", "Gap", "Nodes", "Root Relaxation", "Error_test", "Error_test_Post_Pr"};
  for (int c = 0; c<16; c++){
    columns_name[c] = cl_n[c];
  }
  
  nb_lines = nb_mod*nb_p*nb_dep;
  nb_models= nb_mod;
  nb_part = nb_p;
  nb_depths = nb_dep;
  
  dataset = dts;
  Nmin = Nm;
  alpha = alph;
  time_l = tl;

  content = new string[nb_columns*nb_lines];
}

void statistical_test_matrix::write_line(int p, int d, int model_cpt, string model, solution sol, string error_test, string error_test_pp){
  int cpt = ((p*nb_depths+d-2)*nb_models + model_cpt)*nb_columns;
  
  content[cpt + 0] = dataset;
  content[cpt + 1] = to_string(p);
  content[cpt + 2] = model;
  content[cpt + 3] = to_string(time_l);
  content[cpt + 4] = to_string(d);
  content[cpt + 5] = to_string(alpha);
  content[cpt + 6] = to_string(Nmin);
  
  content[cpt + 7] = to_string(sol.obj);
  content[cpt + 8] = to_string(sol.error_train);
  content[cpt + 9] = to_string(sol.nb_br);
  content[cpt + 10] = to_string(sol.time);
  content[cpt + 11] = to_string(sol.gap);
  content[cpt + 12] = to_string(sol.nodes);
  content[cpt + 13] = to_string(sol.root_rel);
  content[cpt + 14] = error_test;
  content[cpt + 15] = error_test_pp;
}

void statistical_test_matrix::write_result_partition(int p){
  int beg_line = p*nb_depths*nb_models,
    end_line = (p+1)*nb_depths*nb_models;
  string namefile = "../results/"+dataset+"_tl"+to_string(time_l)+"_alp"+to_string(alpha) + "_part" + to_string(p) + "_test.csv";
  write_csv(namefile, beg_line, end_line);
}

learning_test_matrix::learning_test_matrix(int nb_mod, int nb_p, string dts, int D){
  nb_columns = 26;
  columns_name = new string[26];
  string cl_n[26] = { "Dataset", "Partition", "Model", "Time_Limit","DMAX", "Temps","Nb opti","Nb iter",
		      "%err_m1","alpha_m1","best_tree_m1","%err_m2","alpha_m2","best_tree_m2","%err_m3","alpha_m3","best_tree_m3",
		      "%err_m4","alpha_m4","best_tree_m4","%err_m5","alpha_m5","best_tree_m5","%err_m6","alpha_m6","best_tree_m6"};
  for (int c = 0; c<26; c++){
    columns_name[c] = cl_n[c];
  }
  
  nb_lines = nb_mod*nb_p;
  nb_models= nb_mod;
  nb_part = nb_p;
  DMAX = D;
  
  dataset = dts;

  content = new string[nb_columns*nb_lines];
}

void learning_test_matrix::write_line(int p, int model_cpt, string model, double time_l, training_results tr){
  int cpt = (p*nb_models + model_cpt)*nb_columns;
  
  content[cpt + 0] = dataset;
  content[cpt + 1] = to_string(p);
  content[cpt + 2] = model;
  content[cpt + 3] = to_string(time_l);
  content[cpt + 4] = to_string(DMAX);
  content[cpt + 5] = to_string(tr.time);
  content[cpt + 6] = to_string(tr.nb_opti);
  content[cpt + 7] = to_string(tr.nb_iter);
  for (int i = 0; i<6; i++){
    content[cpt + 8 + i*3] = to_string(tr.perc_err[i]);
    content[cpt + 9 + i*3] = to_string(tr.alph[i]);
    content[cpt + 10 + i*3] = tr.best_tree[i];
  }
}

void learning_test_matrix::write_result_partition(int p){
  int beg_line = p*nb_models,
    end_line = (p+1)*nb_models;
  string namefile = "../results/"+dataset+"_DMAX"+to_string(DMAX)+ "_part" + to_string(p) + "_learning_test.csv";
  write_csv(namefile, beg_line, end_line);
}

void test(string dataset_name, int depth_max, int Nmin, float alph, double mu, float time_l){
  GRBEnv env = GRBEnv();
  dataset dt = dataset(dataset_name);
  
  int nb_models = 10,
    nb_part = 10,
    nb_depths = depth_max - 1;

  statistical_test_matrix results = statistical_test_matrix(nb_models, nb_part, nb_depths, dataset_name, Nmin, alph, time_l);
    
  solution sol;
  int model_cpt;

  for (int p = 0; p<nb_part; p++){
    dataset dt_train, dt_test;
    dt.partitionning(dt_train, dt_test);
    //dt.partitionning(dt_train, dt_test,0.2,"partition_"+to_string(p)+".txt");
    
    for (int d = 2; d <= depth_max; d++){
      parameters param = parameters(d,dt_train,alph,mu,0,true,Nmin); // ici on à alpha > 0 et les contraint sum dt <= C , n'existent pas

      model_cpt = 0;

      cout << "OCT\n" ;
      GRBModel oct = GRBModel(env);
      build_model OCT = build_model(oct,dt_train,model_type(baseModel::OCT,true,false,false,false),param);
      sol = OCT.solve(oct,time_l);
      int et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, true);
      int et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"OCT",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "QOCT\n" ;
      GRBModel qoct = GRBModel(env);
      build_model QOCT = build_model(qoct,dt_train,model_type(baseModel::QOCT,true,false,false,false),param);
      sol = QOCT.solve(qoct,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"QOCT",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "FOCT\n" ;
      GRBModel foct = GRBModel(env);
      build_model FOCT = build_model(foct,dt_train,model_type(baseModel::FOCT,true,false,false,false),param);
      sol = FOCT.solve(foct,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"FOCT",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "GOCT\n" ;
      GRBModel goct = GRBModel(env);
      build_model GOCT = build_model(goct,dt_train,model_type(baseModel::GOCT,true,false,false,false),param);
      sol = GOCT.solve(goct,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"GOCT",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "F\n";
      GRBModel f = GRBModel(env);
      build_model F = build_model(f,dt_train,model_type(baseModel::F,true,false,false,false),param);
      sol = F.solve(f,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, false);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"F",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "OCTH\n" ;
      GRBModel octh = GRBModel(env);
      build_model OCTH = build_model(octh,dt_train,model_type(baseModel::OCT,false,false,false,false),param);
      sol = OCTH.solve(octh,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"OCTH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "QOCTH\n" ;
      GRBModel qocth = GRBModel(env);
      build_model QOCTH = build_model(qocth,dt_train,model_type(baseModel::QOCT,false,false,false,false),param);
      sol = QOCTH.solve(qocth,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"QOCTH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "FOCTH\n" ;
      GRBModel focth = GRBModel(env);
      build_model FOCTH = build_model(focth,dt_train,model_type(baseModel::FOCT,false,false,false,false),param);
      sol = FOCTH.solve(focth,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"FOCTH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "GOCTH\n" ;
      GRBModel gocth = GRBModel(env);
      build_model GOCTH = build_model(gocth,dt_train,model_type(baseModel::GOCT,false,false,false,false),param);
      sol = GOCTH.solve(gocth,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"GOCTH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      cout << "FH\n" ;
      GRBModel fh = GRBModel(env);
      build_model FH = build_model(fh,dt_train,model_type(baseModel::F,false,false,false,false),param);
      sol = FH.solve(fh,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, false);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"FH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

    }
    results.write_result_partition(p);
  }
  results.write_csv("../results/"+dataset_name+"_tl"+to_string(time_l)+"_alp"+to_string(alph) + "_test.csv");
}

void learning_test(string dataset_name, int depth_max, double mu){
  GRBEnv env = GRBEnv();
  dataset dt = dataset(dataset_name);

  int Nmin = (int)floor(0.05*dt.I);
  
  int nb_models = 6,
    nb_part = 5;

  learning_test_matrix results = learning_test_matrix(nb_models, nb_part, dataset_name, depth_max);
    
  solution sol;
  int model_cpt;

  time_t now;
  char *date;
  
  for (int p = 0; p<nb_part; p++){
    dataset dt_train, dt_test, dt_validation;
    dt.readPartition(p,dt_train, dt_validation, dt_test);

    model_cpt = 0;

    training_results tr;

    int  Nmin = (int)floor(dt_train.I*0.05);
    double time_limit_univ = 1800.0,
      time_limit_multiv = 300.0;

    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model OCT - start : "<< date<<endl; 
    tr = learning_Bertsimas(dt_train, dt_validation, dt_test, baseModel::OCT, true, depth_max, time_limit_univ, Nmin);
    results.write_line(p,model_cpt,"OCT",time_limit_univ,tr);
    model_cpt += 1;
    
    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model FOCT - start : "<< date<<endl; 
    tr = learning(dt_train, dt_validation, dt_test, baseModel::FOCT, true, depth_max, time_limit_univ, Nmin);
    results.write_line(p,model_cpt,"FOCT",time_limit_univ,tr);
    model_cpt += 1;

    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model F - start : "<< date<<endl; 
    tr = learning(dt_train, dt_validation, dt_test, baseModel::F, true, depth_max, time_limit_univ,0);
    results.write_line(p,model_cpt,"F",time_limit_univ,tr);
    model_cpt += 1;

    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model OCTH - start : "<< date<<endl; 
    tr = learning_Bertsimas(dt_train, dt_validation, dt_test, baseModel::OCT, false, depth_max, time_limit_multiv, Nmin);
    results.write_line(p,model_cpt,"OCTH",time_limit_multiv,tr);
    model_cpt += 1;

    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model FOCTH - start : "<< date<<endl; 
    tr = learning(dt_train, dt_validation, dt_test, baseModel::FOCT, false, depth_max, time_limit_multiv, Nmin);
    results.write_line(p,model_cpt,"FOCTH",time_limit_multiv,tr);
    model_cpt += 1;

    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model FH - start : "<< date<<endl; 
    tr = learning(dt_train, dt_validation, dt_test, baseModel::F, false, depth_max, time_limit_multiv, 0);
    results.write_line(p,model_cpt,"FH",time_limit_multiv,tr);
    model_cpt += 1;

    results.write_result_partition(p);
  }
  
  results.write_csv("../results/"+dataset_name+"_DMAX"+to_string(depth_max) + "_learning_test.csv");
}
 
void testClust(string datasetName){
  GRBEnv env = GRBEnv();
  dataset dt = dataset(datasetName);

  int nbPart = 5;
  
  for (int p=0; p<nbPart; p++){
    dataset dt_train, dt_test, dt_validation;
    dt.readPartition(p,dt_train, dt_validation, dt_test);

    int Nmin = (int)floor(0.05*dt_train.I);

    vector<clustering> clz;
    vector<string> clzNames;
    for (int gm=1; gm<10; gm++){
      clz.push_back(greedyClustering(dt_train,gm/10));
      clzNames.push_back("Greedy"+to_string(gm)+"0%");
    }
    clz.push_back(homogeneousClustering(dt_train, false));
    clzNames.push_back("AlgoLongBary");
    clz.push_back(homogeneousClustering(dt_train, true));
    clzNames.push_back("AlgoLongMedoid");
    clz.push_back(weightedGreedyClustering(dt_train));
    clzNames.push_back("BetterGreedy");
    
    for (int D=2; D<5; D++){
      int Cmin = D,
  	CmaxU = pow(2,D)-1,
	CmaxM = dt.J*(pow(2,D) -1);
      
      model_type mtz[4] = {model_type(baseModel::FOCT,true, false, true, false),
			model_type(baseModel::FOCT,true, false, true, false),
			model_type(baseModel::FOCT,false, false, true, false),
			model_type(baseModel::FOCT,false, false, true, false)};

      parameters paramz[4] = {parameters(D, dt_train, 1/(Cmin*2), 0.00001, Cmin, false, Nmin),
			      parameters(D, dt_train, 1/(CmaxU*2), 0.00001, CmaxU, false, Nmin),
			      parameters(D, dt_train, 1/(Cmin*2), 0.00001, Cmin, false, Nmin),
			      parameters(D, dt_train, 1/(CmaxM*2), 0.00001, CmaxM, false, Nmin)};
      for (int i=0; i<4; i++){
	string typeOfSplit;
	if (mtz[i].univ){
	  typeOfSplit = "UNIV";
	}
	else{
	  typeOfSplit = "MULTIV";
	}
	
	for (int c=0; c<clz.size(); c++){
	  string filename = "../resultsIteratingAlgo/" + dt.name + "_part"+to_string(p)+"_" + clzNames[c]+"_"+ typeOfSplit +"D="+to_string(D) +"_C="+to_string(paramz[i].C)+".txt";
	  solClust solC = decontractingAlgorithm(dt_train, clz[c], mtz[i], paramz[i]);
	  solC.errTr = solC.T.prediction_errors(dt_train);
	  solC.errTst = solC.T.prediction_errors(dt_test);
	  solC.write(filename);
	}
	string filename = "../resultsIteratingAlgo/" + dt.name + "_part"+to_string(p)+"_NoClust_"+ typeOfSplit +"D="+to_string(D) +"_C="+to_string(paramz[i].C)+".txt";

	GRBModel foct = GRBModel(env);
	build_model FOCT = build_model(foct,dt_train,mtz[i],paramz[i]);
	solution sol = FOCT.solve(foct,3600);
     

	fstream file;
	file.open(filename,ios::out);
	file << sol.time <<endl;
	file << sol.T.prediction_errors(dt_train) << endl;
	file << sol.T.prediction_errors(dt_test) << endl;
	file.close();
      }
    }
  }
}
