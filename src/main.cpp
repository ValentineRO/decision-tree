#include "test.h"
#include "testv2.h"
#include <vector>
#include <algorithm>
using namespace std;

int main(){

  /*
  //string datasetNames[4] = {"iris","wine","monk2","haberman"};
  string datasetNames[1] = {"iris"};
  vector<pair<int,int>> DC = {pair<int,int>(2,2),pair<int,int>(2,3),pair<int,int>(3,3),pair<int,int>(3,5),pair<int,int>(3,7),
			      pair<int,int>(4,4),pair<int,int>(4,7),pair<int,int>(4,10)};

  string error = "";
  for (auto dtName: datasetNames){
    cout << dtName << endl;
    dataset dt = dataset(dtName);

    for (int p=0; p<5; p++){
      cout <<"Partition " << p << endl;
      dataset dt_train, dt_validation, dt_test;
      dt.readPartition(p, dt_train, dt_validation, dt_test);

      int  Nmin = (int)floor(dt_train.I*0.05);

      dt_train.computeDists();

      vector<clustering> clz;
      vector<time_t> clzTime;
      time_t t1,t2;
      t1 = time (NULL);
      clz.push_back(hierarchicalClustering(dt_train, 0.25));
      t2 = time (NULL);
      clzTime.push_back(t2-t1);
      t1 = time (NULL);
      clz.push_back(weightedGreedyClustering(dt_train));
      t2 = time (NULL);
      clzTime.push_back(t2-t1);
      t1 = time (NULL);
      clz.push_back(homogeneousClustering(dt_train,true));
      t2 = time (NULL);
      clzTime.push_back(t2-t1);

      for (int d = 0; d<DC.size(); d++){
	cout << "D = " << DC[d].first << " - C = " << DC[d].second << endl;
	for (int s=0; s<2; s++){
	  cout << ((s==0) ? "UNIV" : "MULTIV") << endl;
	  parameters pm1 = parameters(DC[d].first, dt_train, (double)1/(double)(DC[d].second+1), DC[d].second, false, Nmin);
	  model_type mt1 = model_type(baseModel::FOCT, s==0, false, true, false);
	  model_type mt2 = model_type(baseModel::FOCT, s==0, false, true, false, true);

	  string fileName = "./rez/"+dtName+"_part"+to_string(p)+"_D="+to_string(DC[d].first)+"_C="+to_string(DC[d].second);
	  fileName += (s==0) ? "_UNIV" : "_MULTIV";
	  
	  for (int c=0;c<clz.size();c++){
	    cout << "Clustering " << c << endl;
	    dataset dtCl = clz[c].createDt(dt_train);
	    parameters pm2 = parameters(DC[d].first, dtCl, (double)1/(double)(DC[d].second+1), DC[d].second, false, Nmin);

	    cout << "Heuristic solve" << endl;
	    try {
	      solClust solC1 = approxIteratingOTP(dt_train, clz[c], mt1, pm2, false, 1800);
	      solC1.clTime = clzTime[c];
	      solC1.write(fileName+"_cl"+to_string(c)+"_approx.txt");
	    } catch (...){
	      error += fileName+"_cl"+to_string(c)+"_approx.txt \n";
	      
	    }
	   
	    cout << "Exact solve" << endl;
	    try{
	      solClust solC2 = iteratingOTP(dt_train, clz[c], mt2, pm2, false, 1800);
	      solC2.clTime = clzTime[c];
	      solC2.write(fileName+"_cl"+to_string(c)+".txt");
	    } catch (...){
	      error += fileName+"_cl"+to_string(c)+".txt \n";
	    }
	  }
	}
      }
    }
  }

  cout << error;

  */

  /*  
  //string datasetNames[4] = {"iris","wine","monk2","haberman"};
  string datasetNames[1] = {""};
  vector<pair<int,int>> DC = {pair<int,int>(2,2),pair<int,int>(2,3),pair<int,int>(3,3),pair<int,int>(3,5),pair<int,int>(3,7),
			      pair<int,int>(4,4),pair<int,int>(4,7),pair<int,int>(4,10)};

  string error = "";
  for (auto dtName: datasetNames){
    cout << dtName << endl;
    dataset dt = dataset(dtName);

    for (int p=0; p<5; p++){
      cout <<"Partition " << p << endl;
      dataset dt_train, dt_validation, dt_test;
      dt.readPartition(p, dt_train, dt_validation, dt_test);

      int  Nmin = (int)floor(dt_train.I*0.05);

      dt_train.computeDists();

      vector<clustering> clz;
      vector<time_t> clzTime;
      time_t t1,t2;
      t1 = time (NULL);
      clz.push_back(hierarchicalClustering(dt_train, 0.25));
      t2 = time (NULL);
      clzTime.push_back(t2-t1);
      t1 = time (NULL);
      clz.push_back(weightedGreedyClustering(dt_train));
      t2 = time (NULL);
      clzTime.push_back(t2-t1);
      t1 = time (NULL);
      clz.push_back(homogeneousClustering(dt_train,true));
      t2 = time (NULL);
      clzTime.push_back(t2-t1);

      for (int d = 0; d<DC.size(); d++){
	cout << "D = " << DC[d].first << " - C = " << DC[d].second << endl;
	for (int s=0; s<2; s++){
	  cout << ((s==0) ? "UNIV" : "MULTIV") << endl;
	  parameters pm1 = parameters(DC[d].first, dt_train, (double)1/(double)(DC[d].second+1), DC[d].second, false, Nmin);
	  model_type mt1 = model_type(baseModel::FOCT, s==0, false, true, false);
	  model_type mt2 = model_type(baseModel::FOCT, s==0, false, true, false, true);

	  string fileName = "./rez/"+dtName+"_part"+to_string(p)+"_D="+to_string(DC[d].first)+"_C="+to_string(DC[d].second);
	  fileName += (s==0) ? "_UNIV" : "_MULTIV";

	  for (int c=0;c<clz.size();c++){
	    cout << "Clustering " << c << endl;
	    dataset dtCl = clz[c].createDt(dt_train);
	    parameters pm2 = parameters(DC[d].first, dtCl, (double)1/(double)(DC[d].second+1), DC[d].second, false, Nmin);

	    cout << "Heuristic solve" << endl;
	    try{
	      solClust solC1 = approxIteratingOTP(dt_train, clz[c], mt1, pm2, true, 1800);
	      solC1.clTime = clzTime[c];
	      solC1.write(fileName+"_cl"+to_string(c)+"_SansRec_approx.txt");
	    }catch(...){
	      error += fileName+"_cl"+to_string(c)+"_SansRec_approx.txt \n" ;
	    }

	    cout << "Exact solve" << endl;
	    try{
	      solClust solC2 = iteratingOTP(dt_train, clz[c], mt2, pm2, true, 1800);
	      solC2.clTime = clzTime[c];
	      solC2.write(fileName+"_cl"+to_string(c)+"_SansRec.txt");
	    }catch(...){
	      error += fileName+"_cl"+to_string(c)+"_SansRec.txt \n" ;
	    }
	  }
	}
      }
    }
  }

  cout << error;
  */

  
  dataset dt = dataset("wine");
  
  dataset dt_train, dt_validation, dt_test;
  dt.readPartition(0, dt_train, dt_validation, dt_test);

  int  Nmin = (int)floor(dt_train.I*0.05);
  dt_train.computeDists();
  
  clustering cl = hierarchicalClustering(dt_train, 0.25);
//clustering cl = weightedGreedyClustering(dt_train);
  dataset dtCl = cl.createDt(dt_train);

  int D = 2,
    C = 3;
  
  parameters pm = parameters(D, dtCl, (double)1/(double)(C+1), C, false, Nmin);
  model_type mt = model_type(baseModel::FOCT, false, false, true, false, true);

  try {
    //solClust solC = iteratingOTP(dt_train, cl, mt, pm, false, 1800);
    solClust solC = iteratingOTP(dt_train, cl, mt, pm, true, 1800);
    solC.write("youpi.txt");
  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during optimization" << endl;
  }
  
  
  return 0;
}


/*
void getTrainingPerf(){
  string datasetNames[22] = {"ecoli", "iris", "haberman", "monk1", "monk2", "monk3", "seeds", "wine", "balance-scale", "biodeg",
			     "blood-transfusion", "breast-cancer", "car_evaluation", "dermatology", "german", "IndiansDiabetes",
			     "Ionosphere", "tic-tac-toe", "bank_conv", "seismic_bumps", "spambase","Statlog_satellite"};

  string models[6] = {"OCT","FOCT","F","OCTH","FOCTH","FH"};
  
  string filename("../results/learningTrainPerf.csv");
  ofstream file_out;
  file_out.open(filename, ios_base::app);
  file_out << "Dataset,Partition,Profondeur,Model,ErrorTrain" << endl;

  
  for (int dt=0; dt<22; dt++){
    cout << datasetNames[dt] << endl;
    dataset ds = dataset(datasetNames[dt]);
    
    for (int p=0; p<5; p++){
      dataset ds_train, ds_validation, ds_test;
      ds.readPartition(p,ds_train, ds_validation, ds_test);
      
      fstream file;
      file.open("../results/"+datasetNames[dt]+"_DMAX4_part"+to_string(p)+"_learning_test.csv", ios::in);

      string line;
      file >> line; // Dataset,Partition,Model,Time_Limit,Depth,Temps,Nb
      if (line!=""){
	file >> line; // opti,Nb
	file >> line; // iter,%errWithoutPP,alphaWithoutPP,bestTreeWithoutPP,%errWithPP,alphaWithPP,bestTreeWithPP,%errBestOfBoth,alphaBestOfBoth,bestTreeBestOfBoth,
	for (int m=0; m<6; m++){
	  for (int d=2; d<5; d++){
	    file >> line;
	    vector<string> parsedLine = parseLine(line);

	    string treeName;
	    if (m==0 or m==3){
	      treeName = parsedLine[10];
	    }
	    else{
	      treeName = parsedLine[16];
	    }

	    Tree T = Tree(treeName);
	    int trainError = T.prediction_errors(ds_train);
	    file_out << datasetNames[dt] << "," << p << "," << d << "," << models[m] << "," << (float)trainError/(float)ds_train.I << endl;
	    //cout << datasetNames[dt] << "," << p << "," << d << "," << models[m] << "," << (float)trainError/(float)ds_train.I << endl;
	  }
	}
      }
      file.close();
    }
  }
  
  file_out.close();

}

void getLeaves(){
  string datasetNames[22] = {"ecoli", "iris", "haberman", "monk1", "monk2", "monk3", "seeds", "wine", "balance-scale", "biodeg",
			     "blood-transfusion", "breast-cancer", "car_evaluation", "dermatology", "german", "IndiansDiabetes",
			     "Ionosphere", "tic-tac-toe", "bank_conv", "seismic_bumps", "spambase","Statlog_satellite"};

  string models[6] = {"OCT","FOCT","F","OCTH","FOCTH","FH"};
  
  string filename("../results/nbLeaves.csv");
  ofstream file_out;
  file_out.open(filename, ios_base::app);
  file_out << "Dataset,Partition,Profondeur,Model,Leaves" << endl;


  int nbError = 0;
  
  for (int dt=0; dt<22; dt++){
    cout << datasetNames[dt] << endl;
    
    for (int p=0; p<5; p++){
      fstream file;
      try{
	file.open("../results/"+datasetNames[dt]+"_DMAX4_part"+to_string(p)+"_learning_test.csv", ios::in);

	string line;
	file >> line; // Dataset,Partition,Model,Time_Limit,Depth,Temps,Nb
	if (line!=""){
	  file >> line; // opti,Nb
	  file >> line; // iter,%errWithoutPP,alphaWithoutPP,bestTreeWithoutPP,%errWithPP,alphaWithPP,bestTreeWithPP,%errBestOfBoth,alphaBestOfBoth,bestTreeBestOfBoth,
	  for (int m=0; m<6; m++){
	    for (int d=2; d<5; d++){
	      file >> line;
	      vector<string> parsedLine = parseLine(line);

	      string treeName;
	      if (m==0 or m==3){
		treeName = parsedLine[10];
	      }
	      else{
		treeName = parsedLine[16];
	      }

	      Tree T = Tree(treeName);
	      int nbLeaves = 0;
	      for (int l=0; l<T.N+T.L; l++){
		if (T.c[l] != -1){
		  nbLeaves += 1;
		}
	      }

	      file_out << datasetNames[dt] << "," << p << "," << d << "," << models[m] << "," << nbLeaves << endl;
	    }
	  }
	}
	file.close();
      } catch(...){
	nbError += 1;
      }
      
    }
  }
  
  file_out.close();
  cout << nbError<<endl;
}

vector<string> parseLine(string line){
  vector<string> v;
  v.push_back("");
  for (int i=0; i<line.size(); i++){
    if (line[i] == ','){
      v.push_back("");
    }
    else{
      v[v.size()-1] += line[i];
    }
  }
  return v;
}
*/
