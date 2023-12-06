#include "test.h"

extern GurobiEnvironment& gurobiEnv;

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


statistical_test_matrix::statistical_test_matrix(int nb_mod, int nb_p, int nb_dep, string dts, int Nm, float alph, vector<double> tl){
  nb_columns = 17;
  columns_name = new string[17];
  string cl_n[17] = { "Dataset", "Partition", "Model", "Time_Limit", "Profondeur", "Alpha", "Nmin", "Objective", "Error train", "Number of branchings", "Time", "Gap", "Nodes", "Root Relaxation", "Error_test", "Error_test_Post_Pr","ObjectiveEvolution"};
  for (int c = 0; c<17; c++){
    columns_name[c] = cl_n[c];
  }
  
  nb_lines = nb_mod*nb_p*nb_dep;
  nb_models= nb_mod;
  nb_part = nb_p;
  nb_depths = nb_dep;
  
  dataset = dts;
  Nmin = Nm;
  alpha = alph;
  timeL = tl[0] + tl[1] + tl[2] + tl[3] + tl[4] + tl[5]; // ya que 6 temps différents

  content = new string[nb_columns*nb_lines];
}

void statistical_test_matrix::write_line(int p, int d, int model_cpt, string model, solution sol, string error_test, string error_test_pp){
  int cpt = ((p*nb_depths+d-2)*nb_models + model_cpt)*nb_columns;
  
  content[cpt + 0] = dataset;
  content[cpt + 1] = to_string(p);
  content[cpt + 2] = model;
  content[cpt + 3] = to_string(timeL);
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
  //content[cpt + 16] = sol.objEvo;
}

void statistical_test_matrix::write_result_partition(int p){
  int beg_line = p*nb_depths*nb_models,
    end_line = (p+1)*nb_depths*nb_models;
  string namefile = "../results/"+dataset+"_tl"+to_string(timeL)+"_alp"+to_string(alpha) + "_part" + to_string(p) + "_test.csv";
  write_csv(namefile, beg_line, end_line);
}

learning_test_matrix::learning_test_matrix(int nb_mod, int nb_p, string dts, int nb_d){
  nb_columns = 17;
  columns_name = new string[17];
  string cl_n[17] = { "Dataset", "Partition", "Model", "Time_Limit","Depth", "Temps","Nb opti","Nb iter",
		      "%errWithoutPP","alphaWithoutPP","bestTreeWithoutPP","%errWithPP","alphaWithPP","bestTreeWithPP",
		      "%errBestOfBoth","alphaBestOfBoth","bestTreeBestOfBoth"};
  for (int c = 0; c<17; c++){
    columns_name[c] = cl_n[c];
  }
  
  nb_lines = nb_mod*nb_p*nb_d;
  nb_models= nb_mod;
  nb_part = nb_p;
  nb_depths = nb_d;
  
  dataset = dts;

  content = new string[nb_columns*nb_lines];
}

void learning_test_matrix::write_lines(int p, int model_cpt, string model, bool thereIsPP, double time_l, training_results tr){
  int cpt = (p*nb_models*nb_depths + model_cpt*nb_depths)*nb_columns;

  for (int d=0; d<nb_depths; d++){
    content[cpt + nb_columns*d + 0] = dataset;
    content[cpt + nb_columns*d + 1] = to_string(p);
    content[cpt + nb_columns*d + 2] = model;
    content[cpt + nb_columns*d + 3] = to_string(time_l);
    content[cpt + nb_columns*d + 4] = to_string(d+2);
    content[cpt + nb_columns*d + 5] = to_string(tr.time[d]);
    content[cpt + nb_columns*d + 6] = to_string(tr.nbOpti[d]);
    content[cpt + nb_columns*d + 7] = to_string(tr.nbIter[d]);
   
    content[cpt + nb_columns*d + 8] = to_string(tr.percErrWithoutPP[d]);
    content[cpt + nb_columns*d + 9] = to_string(tr.alphWithoutPP[d]);
    content[cpt + nb_columns*d + 10] = tr.bestTreeWithoutPP[d];

    if (thereIsPP){
      content[cpt + nb_columns*d + 11] = to_string(tr.percErrWithPP[d]);
      content[cpt + nb_columns*d + 12] = to_string(tr.alphWithPP[d]);
      content[cpt + nb_columns*d + 13] = tr.bestTreeWithPP[d];

      content[cpt + nb_columns*d + 14] = to_string(tr.percErrBestOfBoth[d]);
      content[cpt + nb_columns*d + 15] = to_string(tr.alphBestOfBoth[d]);
      content[cpt + nb_columns*d + 16] = tr.bestTreeBestOfBoth[d];
    }
  }
}

void learning_test_matrix::write_result_partition(int p){
  int beg_line = p*nb_models*nb_depths,
    end_line = (p+1)*nb_models*nb_depths;
  string namefile = "../results/"+dataset+"_DMAX"+to_string(nb_depths+1)+ "_part" + to_string(p) + "_learning_test.csv";
  write_csv(namefile, beg_line, end_line);
}

/*
void test(string dataset_name){
  GRBEnv env = GRBEnv();
  dataset dt = dataset(dataset_name);

  int depth_max = 4;
  int Nmin = 0;
  float alph = 0.1;

  vector<double> time_l = {60,60,180,300,600,600};  // ca veut dire au bout de {60,120,300,600,1200,1800} secondes
  
  int nb_models = 10,
    nb_part = 5,
    nb_depths = depth_max - 1;

  statistical_test_matrix results = statistical_test_matrix(nb_models, nb_part, nb_depths, dataset_name, Nmin, alph, time_l);
    
  solution sol;
  int model_cpt;

  time_t now;
  char* date;

  for (int p = 0; p<nb_part; p++){
    dataset dt_train, dt_validation, dt_test;
    dt.readPartition(p,dt_train, dt_validation, dt_test);
    
    for (int d = 2; d <= depth_max; d++){
      parameters param = parameters(d,dt_train,alph,-1,true,Nmin); // ici on à alpha > 0 et les contraint sum dt <= C , n'existent pas

      model_cpt = 0;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " OCT - start : "<< date<<endl;
      GRBModel oct = GRBModel(env);
      build_model OCT = build_model(oct,dt_train,model_type(baseModel::OCT,true,false,false,false),param);
      sol = OCT.solve(oct,0.0,time_l);
      int et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, true);
      int et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"OCT",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " QOCT - start : "<< date<<endl;
      GRBModel qoct = GRBModel(env);
      build_model QOCT = build_model(qoct,dt_train,model_type(baseModel::QOCT,true,false,false,false),param);
      sol = QOCT.solve(qoct,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"QOCT",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " FOCT - start : "<< date<<endl;
      GRBModel foct = GRBModel(env);
      build_model FOCT = build_model(foct,dt_train,model_type(baseModel::FOCT,true,false,false,false),param);
      sol = FOCT.solve(foct,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"FOCT",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " GOCT - start : "<< date<<endl;
      GRBModel goct = GRBModel(env);
      build_model GOCT = build_model(goct,dt_train,model_type(baseModel::GOCT,true,false,false,false),param);
      sol = GOCT.solve(goct,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"GOCT",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " F - start : "<< date<<endl;
      GRBModel f = GRBModel(env);
      build_model F = build_model(f,dt_train,model_type(baseModel::F,true,false,false,false),param);
      sol = F.solve(f,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_b(dt_train, false);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"F",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " OCTH - start : "<< date<<endl;
      GRBModel octh = GRBModel(env);
      build_model OCTH = build_model(octh,dt_train,model_type(baseModel::OCT,false,false,false,false),param);
      sol = OCTH.solve(octh,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"OCTH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " QOCTH - start : "<< date<<endl;
      GRBModel qocth = GRBModel(env);
      build_model QOCTH = build_model(qocth,dt_train,model_type(baseModel::QOCT,false,false,false,false),param);
      sol = QOCTH.solve(qocth,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"QOCTH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " FOCTH - start : "<< date<<endl;
      GRBModel focth = GRBModel(env);
      build_model FOCTH = build_model(focth,dt_train,model_type(baseModel::FOCT,false,false,false,false),param);
      sol = FOCTH.solve(focth,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"FOCTH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " GOCTH - start : "<< date<<endl;
      GRBModel gocth = GRBModel(env);
      build_model GOCTH = build_model(gocth,dt_train,model_type(baseModel::GOCT,false,false,false,false),param);
      sol = GOCTH.solve(gocth,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, true);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"GOCTH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " FH - start : "<< date<<endl;
      GRBModel fh = GRBModel(env);
      build_model FH = build_model(fh,dt_train,model_type(baseModel::F,false,false,false,false),param);
      sol = FH.solve(fh,0.0,time_l);
      et = sol.T.prediction_errors(dt_test);
      sol.T.post_processing_a_b(dt_train, false);
      et_pp = sol.T.prediction_errors(dt_test);
      results.write_line(p,d,model_cpt,"FH",sol,to_string(et),to_string(et_pp));
      model_cpt += 1;

    }
    results.write_result_partition(p);
  }
  results.write_csv("../results/"+dataset_name+"_tl"+to_string(results.timeL)+"_alp"+to_string(alph) + "_test.csv");
}
*/

void learning_test(string dataset_name, int depth_max, bool doOCT){
  //GRBEnv env = GRBEnv();
  GRBEnv& env = gurobiEnv.getEnvironment();
  dataset dt = dataset(dataset_name);

  int nb_models = 6,
    nb_part = 5,
    nb_depth = depth_max -1;

  if (!doOCT){
    nb_models = 4;
  }

  learning_test_matrix results = learning_test_matrix(nb_models, nb_part, dataset_name, nb_depth);
    
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

    if (doOCT){
      now = time(0);
      date = ctime(& now);
      cout << endl;
      cout << "Partition "<< p << " / " << nb_part << " - Model OCT - start : "<< date<<endl; 
      tr = learning_Bertsimas(dt_train, dt_validation, dt_test, baseModel::OCT, true, depth_max, time_limit_univ, Nmin);
      results.write_lines(p,model_cpt,"OCT",false,time_limit_univ,tr);
      model_cpt += 1;
    }
    
    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model FOCT - start : "<< date<<endl; 
    tr = learning(dt_train, dt_validation, dt_test, baseModel::FOCT, true, depth_max, time_limit_univ, Nmin);
    results.write_lines(p,model_cpt,"FOCT",true,time_limit_univ,tr);
    model_cpt += 1;

    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model F - start : "<< date<<endl; 
    tr = learning(dt_train, dt_validation, dt_test, baseModel::F, true, depth_max, time_limit_univ,0);
    results.write_lines(p,model_cpt,"F",true,time_limit_univ,tr);
    model_cpt += 1;

    
    if (doOCT){
      now = time(0);
      date = ctime(& now);
      cout << endl;
      cout << "Partition "<< p << " / " << nb_part << " - Model OCTH - start : "<< date<<endl; 
      tr = learning_Bertsimas(dt_train, dt_validation, dt_test, baseModel::OCT, false, depth_max, time_limit_multiv, Nmin);
      results.write_lines(p,model_cpt,"OCTH",false,time_limit_multiv,tr);
      model_cpt += 1;
    }

    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model FOCTH - start : "<< date<<endl; 
    tr = learning(dt_train, dt_validation, dt_test, baseModel::FOCT, false, depth_max, time_limit_multiv, Nmin);
    results.write_lines(p,model_cpt,"FOCTH",true,time_limit_multiv,tr);
    model_cpt += 1;

    now = time(0);
    date = ctime(& now);
    cout << endl;
    cout << "Partition "<< p << " / " << nb_part << " - Model FH - start : "<< date<<endl; 
    tr = learning(dt_train, dt_validation, dt_test, baseModel::F, false, depth_max, time_limit_multiv, 0);
    results.write_lines(p,model_cpt,"FH",true,time_limit_multiv,tr);
    model_cpt += 1;

    results.write_result_partition(p);
  }
  results.write_csv("../results/"+dataset_name+"_DMAX"+to_string(depth_max) + "_learning_test.csv");
}

/*
void testClustAncien(string datasetName, float time_l){
  GRBEnv env = GRBEnv();
  dataset dt = dataset(datasetName);

  int nbPart = 5;
  int nbClustering = 12;

  string column_names[8] = {"numPart","clustType","time","% red","homogeneity","exclusion","consistency","meanDist"};

  info_matrix iM = info_matrix(8,nbPart*nbClustering,column_names);
  
  for (int p=0; p<nbPart; p++){
    cout << "################# Partition " << p << endl;
    dataset dt_train, dt_test, dt_validation;
    dt.readPartition(p,dt_train, dt_validation, dt_test);

    dt_train.computeDists();

    int Nmin = (int)floor(0.05*dt_train.I);

    vector<clustering> clz;
    vector<time_t> clzTimes;
    vector<string> clzNames;

    for (int gm=1; gm<10; gm++){
      time_t t1 = time (NULL);
      clz.push_back(hierarchicalClustering(dt_train,(float)gm/10));
      time_t t2 = time (NULL);
      clzTimes.push_back(t2-t1);
      clzNames.push_back("Greedy"+to_string(gm)+"0%");
    }
    time_t t1 = time (NULL);
    clz.push_back(homogeneousClustering(dt_train, false));
    time_t t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("AlgoLongBary");

    t1 = time (NULL);
    clz.push_back(homogeneousClustering(dt_train, true));
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("AlgoLongMedoid");
    

    t1 = time (NULL);
    clz.push_back(weightedGreedyClustering(dt_train));
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("BetterGreedy");

    for (int c=0; c<clz.size(); c++){
      string line[8];
      line[0] = to_string(p);
      line[1] = clzNames[c];
      line[2] = to_string(clzTimes[c]);
      line[3] = to_string((float)clz[c].clusters.size()/dt_train.I);
      line[4] = to_string(clz[c].computeHomogeneity(dt_train));
      line[5] = to_string(clz[c].computeExclusion(dt_train));
      line[6] = to_string(clz[c].computeConsistency(dt_train));
      line[7] = to_string(clz[c].computeMeanDist(dt_train));

      iM.write_line(p*nbClustering + c, line);
    }

    for (int D=2; D<5; D++){
      int Cmin = D,
  	CmaxU = pow(2,D)-1,
	CmaxM = dt.J*(pow(2,D) -1);
      
      model_type mtz[4] = {model_type(baseModel::FOCT,true, false, true, false),
			model_type(baseModel::FOCT,true, false, true, false),
			model_type(baseModel::FOCT,false, false, true, false),
			model_type(baseModel::FOCT,false, false, true, false)};

      parameters paramz[4] = {parameters(D, dt_train, (double)1/(Cmin+1), Cmin, false, Nmin),
			      parameters(D, dt_train, (double)1/(CmaxU+1), CmaxU, false, Nmin),
			      parameters(D, dt_train, (double)1/(Cmin+1), Cmin, false, Nmin),
			      parameters(D, dt_train, (double)1/(CmaxM+1), CmaxM, false, Nmin)};
      for (int i=0; i<4; i++){
	string typeOfSplit;
	if (mtz[i].univ){
	  typeOfSplit = "UNIV";
	}
	else{
	  typeOfSplit = "MULTIV";
	}
	
	for (int c=0; c<clz.size(); c++){
	  cout << dt.name + "_part"+to_string(p)+"_" + clzNames[c]+"_"+ typeOfSplit +"_D="+to_string(D) +"_C="+to_string(paramz[i].C) <<endl;
	  string filename = "../resultsIteratingAlgo/" + dt.name + "/"+ dt.name + "_part"+to_string(p)+"_" + clzNames[c]+"_"+ typeOfSplit +"_D="+to_string(D) +"_C="+to_string(paramz[i].C)+".txt";
	  solClust solC = approxIteratingOTP(dt_train, clz[c], mtz[i], paramz[i],time_l);
	  solC.addClusteringTime(clzTimes[c]);
	  solC.errTr = solC.finalTree.prediction_errors(dt_train);
	  solC.errTst = solC.finalTree.prediction_errors(dt_test);
	  solC.write(filename);
	}
	string filename = "../resultsIteratingAlgo/" + dt.name +"/"+ dt.name+ "_part"+to_string(p)+"_NoClust_"+ typeOfSplit +"_D="+to_string(D) +"_C="+to_string(paramz[i].C)+".txt";

	GRBModel foct = GRBModel(env);
	t1 = time (NULL);
	build_model FOCT = build_model(foct,dt_train,mtz[i],paramz[i]);
	solution sol = FOCT.solve(foct,time_l);
	t2 = time (NULL);
     
	fstream file;
	file.open(filename,ios::out | std::ofstream::trunc);
	file << t2-t1 <<endl;
	file << sol.T.prediction_errors(dt_train) << endl;
	file << sol.T.prediction_errors(dt_test) << endl;
	file.close();
      }
    }
  }
  iM.write_csv("../resultsIteratingAlgo/"+datasetName+"/clusteringStats.csv");
}

void testClust(string datasetName, float time_l, bool useAlgoLong){
  GRBEnv env = GRBEnv();
  dataset dt = dataset(datasetName);

  int nbPart = 5;
  int nbClustering = 11+(int)useAlgoLong*2;

  string column_names[8] = {"numPart","clustType","time","% red","homogeneity","exclusion","consistency","meanDist"};

  info_matrix iM = info_matrix(8,nbPart*nbClustering,column_names);
  
  for (int p=0; p<nbPart; p++){
    cout << "################# Partition " << p << endl;
    dataset dt_train, dt_test, dt_validation;
    dt.readPartition(p,dt_train, dt_validation, dt_test);

    dt_train.computeDists();

    int Nmin = (int)floor(0.05*dt_train.I);

    vector<clustering> clz;
    vector<time_t> clzTimes;
    vector<string> clzNames;


    time_t t1, t2;
    for (int i=1; i<10; i++){
      t1 = time (NULL);
      clz.push_back(hierarchicalClustering(dt_train,(float)i*0.1));
      t2 = time (NULL);
      clzTimes.push_back(t2-t1);
      clzNames.push_back("Greedy"+to_string(i)+"0%");
    }

    if (useAlgoLong){
      t1 = time (NULL);
      clz.push_back(homogeneousClustering(dt_train, false));
      t2 = time (NULL);
      clzTimes.push_back(t2-t1);
      clzNames.push_back("AlgoLongBary");

      t1 = time (NULL);
      clz.push_back(homogeneousClustering(dt_train, true));
      t2 = time (NULL);
      clzTimes.push_back(t2-t1);
      clzNames.push_back("AlgoLongMedoid");
    }

    t1 = time (NULL);
    clz.push_back(weightedGreedyClustering(dt_train));
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("BetterGreedy");

    t1 = time (NULL);
    clz.push_back(hierarchicalClustering(dt_train,(float)clz[nbClustering-2].clusters.size()/dt_train.I));
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("GreedyCommeBetterGreedy");

    for (int c=0; c<clz.size(); c++){
      string line[8];
      line[0] = to_string(p);
      line[1] = clzNames[c];
      line[2] = to_string(clzTimes[c]);
      line[3] = to_string((float)clz[c].clusters.size()/dt_train.I);
      line[4] = to_string(clz[c].computeHomogeneity(dt_train));
      line[5] = to_string(clz[c].computeExclusion(dt_train));
      line[6] = to_string(clz[c].computeConsistency(dt_train));
      line[7] = to_string(clz[c].computeMeanDist(dt_train));

      iM.write_line(p*nbClustering + c, line);
    }

    for (int D=3; D<5; D++){
      int Cmin = D,
  	CmaxU = pow(2,D)-1,
	CmaxM = dt.J*(pow(2,D) -1);
      
      model_type mtz[4] = {model_type(baseModel::FOCT,true, false, true, false),
			model_type(baseModel::FOCT,true, false, true, false),
			model_type(baseModel::FOCT,false, false, true, false),
			model_type(baseModel::FOCT,false, false, true, false)};

      parameters paramz[4] = {parameters(D, dt_train, (double)1/(Cmin+1), Cmin, false, Nmin),
			      parameters(D, dt_train, (double)1/(CmaxU+1), CmaxU, false, Nmin),
			      parameters(D, dt_train, (double)1/(Cmin+1), Cmin, false, Nmin),
			      parameters(D, dt_train, (double)1/(CmaxM+1), CmaxM, false, Nmin)};
      for (int i=0; i<4; i++){
	string typeOfSplit;
	if (mtz[i].univ){
	  typeOfSplit = "UNIV";
	}
	else{
	  typeOfSplit = "MULTIV";
	}
	
	for (int c=0; c<clz.size(); c++){
	  cout << dt.name + "_part"+to_string(p)+"_" + clzNames[c]+"_"+ typeOfSplit +"_D="+to_string(D) +"_C="+to_string(paramz[i].C) <<endl;
	  string filename = "../resultsIteratingAlgo/" + dt.name + "/"+ dt.name + "_part"+to_string(p)+"_" + clzNames[c]+"_"+ typeOfSplit +"_D="+to_string(D) +"_C="+to_string(paramz[i].C)+".txt";
	  solClust solC = approxIteratingOTP(dt_train, clz[c], mtz[i], paramz[i],time_l);
	  solC.addClusteringTime(clzTimes[c]);
	  solC.errTr = solC.finalTree.prediction_errors(dt_train);
	  solC.errTst = solC.finalTree.prediction_errors(dt_test);
	  solC.write(filename);
	}
	string filename = "../resultsIteratingAlgo/" + dt.name +"/"+ dt.name+ "_part"+to_string(p)+"_NoClust_"+ typeOfSplit +"_D="+to_string(D) +"_C="+to_string(paramz[i].C)+".txt";

	CART CARTalg = CART(dt_train, paramz[i].D, paramz[i].C, paramz[i].Nmin);
	Tree CARTtree = CARTalg.getTree(dt_train);

	cout << dt.name + "_part"+to_string(p)+"_NoClust_"+ typeOfSplit +"_D="+to_string(D) +"_C="+to_string(paramz[i].C) << endl;
	GRBModel foct = GRBModel(env);
	t1 = time (NULL);
	build_model FOCT = build_model(foct,dt_train,mtz[i],paramz[i]);
	FOCT.add_warmstart(CARTtree, dt_train);
	solution sol = FOCT.solve(foct,time_l);
	t2 = time (NULL);
     
	fstream file;
	file.open(filename,ios::out | std::ofstream::trunc);
	file << t2-t1 <<endl;
	file << sol.T.prediction_errors(dt_train) << endl;
	file << sol.T.prediction_errors(dt_test) << endl;
	file.close();
      }
    }
  }
  iM.write_csv("../resultsIteratingAlgo/"+datasetName+"/clusteringStats.csv");
}
*/

void createClusterings(string datasetName){
  dataset dt = dataset(datasetName);

  int nbPart = 5,
    nbClustering;

  vector<string> typeOfClustering;

  if (dt.I/2<1000){
    nbClustering = 7;
    typeOfClustering = {"hierarchicalClustering5", "hierarchicalClustering10","hierarchicalClustering25",
			"betterGreedyClustering", "hierarchicalClusteringGprop",
			"barycenterClustering", "medoidClustering" };
  }
  else{
    nbClustering = 6;
    typeOfClustering = {"hierarchicalClustering1","hierarchicalClustering5", "hierarchicalClustering10",
			"hierarchicalClustering25","betterGreedyClustering", "hierarchicalClusteringGprop"};
  }
  //string column_names[8] = {"Partition","clustType","time","%red","homogeneity","exclusion","consistency","meanDist"};
  string column_names[8] = {"Partition","clustType","consistency"};
  
  //info_matrix iM = info_matrix(8,nbPart*nbClustering,column_names);
  info_matrix iM = info_matrix(3,nbPart*nbClustering,column_names);
  
  for (int p=0; p<nbPart; p++){
    dataset dt_train, dt_test, dt_validation;
    dt.readPartition(p,dt_train, dt_validation, dt_test);

    dt_train.computeDists();

    vector<clustering> clz;
    vector<time_t> clzTimes;
    vector<string> clzNames;

    time_t t1, t2;

    if (dt.I/2>=1000){
      t1 = time (NULL);
      clz.push_back(pythonHierarchicalClustering(dt_train, 0.01));
      t2 = time (NULL);
      clzTimes.push_back(t2-t1);
      clzNames.push_back("hierarchicalClustering1");
    }

    t1 = time (NULL);
    if (dt.I<300){
      clz.push_back(hierarchicalClustering(dt_train, 0.05));
    }
    else{
      clz.push_back(pythonHierarchicalClustering(dt_train, 0.05));
    }
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("hierarchicalClustering5");

    t1 = time (NULL);
    clz.push_back(pythonHierarchicalClustering(dt_train, 0.10));
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("hierarchicalClustering10");

    t1 = time (NULL);
    clz.push_back(pythonHierarchicalClustering(dt_train, 0.25));
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("hierarchicalClustering25");

    t1 = time (NULL);
    clz.push_back(weightedGreedyClustering(dt_train));
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("betterGreedyClustering");

    float percRed = (float)clz.back().clusters.size()/(float)dt_train.I;
    
    t1 = time (NULL);
    clz.push_back(pythonHierarchicalClustering(dt_train, percRed));
    t2 = time (NULL);
    clzTimes.push_back(t2-t1);
    clzNames.push_back("hierarchicalClusteringGprop");

    if (dt.I/2<1000){
      t1 = time (NULL);
      clz.push_back(homogeneousClustering(dt_train));
      t2 = time (NULL);
      clzTimes.push_back(t2-t1);
      clzNames.push_back("barycenterClustering");

      t1 = time (NULL);
      clz.push_back(homogeneousClustering(dt_train,true));
      t2 = time (NULL);
      clzTimes.push_back(t2-t1);
      clzNames.push_back("medoidClustering");
    }


    for (int c=0; c<clz.size(); c++){
      string line[3];
      line[0] = to_string(p);
      line[1] = clzNames[c];
      line[2] = to_string(clz[c].computeConsistency(dt_train));
      //line[3] = to_string((float)clz[c].clusters.size()/dt_train.I);
      //line[4] = to_string(clz[c].computeHomogeneity(dt_train));
      //line[5] = to_string(clz[c].computeExclusion(dt_train));
      //line[6] = to_string(clz[c].computeConsistency(dt_train));
      //line[7] = to_string(clz[c].computeMeanDist(dt_train));

      iM.write_line(p*nbClustering + c, line);

      /*
      clz[c].write("../TreesAndPartitions/" + datasetName + "/part" + to_string(p) + "_" + clzNames[c] + "_NUMBERS.txt");
      dataset clDt = clz[c].createDt(dt_train, c == 3);
      clDt.writeDataset("../TreesAndPartitions/" + datasetName + "/part" + to_string(p) + "_" + clzNames[c] + ".txt");
      */
    }
  }
  //iM.write_csv("../results/clusteringStats_"+ datasetName +".csv");
  iM.write_csv("../results/consistencyStats_"+ datasetName +".csv");
}
