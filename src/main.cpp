#include "test.h"
#include <vector>
#include <algorithm>
using namespace std;

int main(){
  
  //string dt_name[5] = {"iris","wine","dermatology","blood_donation","breast_cancer"};

  /*
  try {
    testClust("iris",300);
  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  */

  dataset dt = dataset("iris");
  dataset dt_train, dt_test, dt_validation;
  dt.readPartition(4,dt_train, dt_validation, dt_test);

  dt_train.computeDists();

  int Nmin = (int)floor(0.05*dt_train.I);

  model_type mt = model_type(baseModel::FOCT,false, false, true, false);
  parameters param = parameters(3, dt_train, (double)1/(3+1), 3, false, Nmin);

  GRBEnv env = GRBEnv();
  GRBModel foct = GRBModel(env);
  build_model FOCT = build_model(foct,dt_train,mt,param);
  solution sol = FOCT.solve(foct,300);
  sol.T.write_tree("uuuuuuuuuue.txt");
  cout << sol.T.prediction_errors(dt_train) << endl;

  clustering cl = homogeneousClustering(dt_train);
  solClust solC = iteratingOTP(dt_train, cl, mt, param,300);

  solC.T.write_tree("ooooooooooooooue.txt");
  cout << solC.T.prediction_errors(dt_train) << endl;

  /* TEST de l'utilisation de la PL pour la création de clusters
  dataset dt = dataset("I100K3P10_inst2");
  dt.computeDists();

  clustering cl1 = homogeneousClustering(dt, false);
  clustering cl2 = optimalClustering(dt, 50);

  cl1.showClusters();
  cl2.showClusters();
  */

  /*
  dataset dt = dataset("iris");
  dataset dt_train, dt_validation, dt_test;
  dt.readPartition(0,dt_train, dt_validation, dt_test);
  dt_train.computeDists();
  
  int Nmin = (int)floor(0.05*dt_train.I);

  model_type mt = model_type(baseModel::FOCT, false, false, true, false);
  parameters param = parameters(4, dt_train, (double)1/(4*2), 4, false, Nmin);

  GRBEnv env = GRBEnv();
  GRBModel foct = GRBModel(env);
  build_model FOCT = build_model(foct,dt_train,mt,param);
  solution sol = FOCT.solve(foct,300);
  sol.T.write_tree("uuuuuuuuuue.txt");

  cout << sol.T.prediction_errors(dt_train) << endl;

  int* z = new int[dt_train.I*sol.T.L];
  FOCT.get_z(z);

  int* leaves = new int[dt_train.I];
  sol.T.predict_leaves(dt_train, leaves);

  for (int i=0; i< dt_train.I; i++){
    if (z[i*sol.T.L+leaves[i]] != 1){
      int correct_leaf;
      for (int l=0; l<sol.T.L; l++){
	if (z[i*sol.T.L+l] == 1){
	  correct_leaf = l;
	}
      }
      cout << i << " is supposed to go in " << correct_leaf << " but goes in " << leaves[i] << endl;
      for (int j=0; j<dt.J; j++){
	cout << "x" << j << " : " << dt_train.X[i*dt.J+j] << endl;
      }
    }
  }
 
  clustering cl = hierarchicalClustering(dt_train, 0.4);
  solClust solC = iteratingOTP(dt_train, cl, mt, param,300);
  */
  

  return 0;
}
