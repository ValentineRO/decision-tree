#include "test.h"
#include "testv2.h"
#include <vector>
#include <algorithm>
using namespace std;

int main(){
  dataset dt = dataset("test");
  clustering cl = clustering(dt, "clTest.txt", "bonsoir");

  dataset dtCl = cl.createDt(dt);

  model_type mt = model_type(baseModel::FOCT, true, false, true, false, true);
  parameters p = parameters(2, dtCl, 0.1, 3, false, 1);

  GRBEnv env = GRBEnv();
  GRBModel m = GRBModel(env);

  build_model md = build_model(m, dtCl, cl, mt, p);

  solution sol = md.solve2(m);

  
  /*
  string datasetName = "wine";
  learningWithClustering(datasetName);
  */

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
