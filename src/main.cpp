#include "test.h"
#include "testv2.h"
#include <vector>
#include <algorithm>
using namespace std;

int main(){

  string dataset_name = "iris";
  learningWithClustering(dataset_name);

  /*  
  try {
    string datasetName = "iris";
    int p=0;
    dataset dt = dataset(datasetName);
    dataset dt_train, dt_test, dt_validation;
    dt.readPartition(p,dt_train, dt_validation, dt_test);

    int Nmin = (int)floor(0.05*dt_train.I),
      D = 3,
      C = 3;
  
    dt_train.computeDists();
    clustering cl = hierarchicalClustering(dt_train, 0.1);
    //clustering cl = weightedGreedyClustering(dt_train);
    //cl.showReduction(dt_train.I);
  
    parameters param = parameters(D,dt_train,(1/(double)(C*2)),C,false,Nmin);

    //CART ogCART = CART(dt, D, C, Nmin);
    //Tree OGcartTree = ogCART.getTree(dt_train);

    //OGcartTree.write_tree("iteration0.txt");
  
    model_type modelt = model_type(baseModel::FOCT, true, false, true, false);

    solClust solC = approxIteratingOTP(dt_train, cl, modelt, param, 1800, Tree());
    //solClust solC = iteratingCART(dt_train, cl, param);
  }catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
    }*/
    
  /*
  for (int p=0; p<5; p++){
    cout << "########################### Partition " << p << endl;
    
    if (p==0){freopen("bsTestA.txt", "a", stdout);}
    if (p==1){freopen("bsTestB.txt", "a", stdout);}
    if (p==2){freopen("bsTestC.txt", "a", stdout);}
    if (p==3){freopen("bsTestD.txt", "a", stdout);}
    if (p==4){freopen("bsTestE.txt", "a", stdout);}
    
 
    dataset dt_train, dt_test, dt_validation;
    dt.readPartition(p,dt_train, dt_validation, dt_test);

    int Nmin = (int)floor(0.05*dt_train.I),
      D = 4,
      C = 4;

    string longClName ="part" + to_string(p) + "_betterGreedyClustering_NUMBERS";
    clustering cl = clustering(dt_train, longClName, "betterGreedyClustering");

    dataset dtCl = cl.createDt(dt_train);

    parameters param = parameters(D,dtCl,0,C,false,Nmin);
    model_type modelt = model_type(baseModel::FOCT, true, false, true, false);

    GRBModel md = GRBModel(env);
    build_model model = build_model(md, dtCl, modelt, param);
    string CARTnamefile = "../TreesAndPartitions/balance-scale/CART_part"+to_string(dt_train.partition)+"_" + cl.name + "_D"+to_string(D)+"_C"+to_string(C)  + "_Nmin.txt";
    Tree warmstartTree = Tree(CARTnamefile);
    model.add_warmstart(warmstartTree, dtCl);

    
    float t = 0;
    int sol_lim = 2;
    solution sol;
    while (t<150){
      md.set("SolutionLimit",to_string(sol_lim));
      sol = model.solve(md, 300-t);
      t += sol.time;
      sol_lim += 1;
    }
    //freopen("/dev/tty", "w", stdout);
    }*/

  /*

  try {
    string dataset_name = "balance-scale";
    createClusterings(dataset_name, true);
    
  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
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
