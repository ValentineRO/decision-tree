#include "data_formatting.h"

void dataset::compute_L_hat(){
  int* count = new int[K];
  int m = Y[0];
  for (int k=0; k<K; k++){
    count[k] = 0;
  }
  for (int i=0; i<I; i++){
    count[Y[i]] += 1;
    if (count[Y[i]] > count[m]){
      m = Y[i];
    }
  }
  L_h = I - count[m];  
}

dataset::dataset(){
  X = {};
  I = 0;
  J = 0;
  K = 0;

  weightedPoints = 0;
}

dataset::dataset(int I_, int J_, int K_) {
  I = I_;
  J = J_;
  K = K_;
  X = new double[I*J];
  Y = new int[I];

  weightedPoints = 0;
}

dataset::dataset(string namefile) {
  fstream file;
  name = namefile;
  partition = -1;
  file.open("../data/"+namefile+".txt", ios::in); // ios::out pour écrire et ios::in pour lire
  file >> I;
  file >> J;
  file >> K;
  //cout << "I : " << I << " J : " << J << " K: " << K << "\n";
  X = new double[I * J];
  int k, l;
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      file >> X[i * J + j];
      //cout << "xij : " << X[i * J + j] << "\n";
    }
  }
  Y = new int[I];
  for (int i = 0; i < I; i++) {
    file >> Y[i];
    //cout << "y : " << Y[i] << "\n";
  }
  file.close();
  compute_mu();
  compute_L_hat();

  weightedPoints = 0;
}

void dataset::compute_mu(){
  mu_vect = new double[J];
  mu_min = 1;
  mu_max = 0;

  for (int j=0; j<J; j++){
    mu_vect[j] = 1;
    for (int i1=0; i1<I; i1++){
      for (int i2=i1+1; i2<I; i2++){
	double space = abs(X[i1*J+j]-X[i2*J+j]);
	if (space != 0.0){
	  mu_vect[j] = min(mu_vect[j],space);
	}
      }
    }
    mu_min = min(mu_min,mu_vect[j]);
    mu_max = max(mu_max, mu_vect[j]);
  }
}

void dataset::computeDists(){
  dists = new double[I*I];

  for (int i1=0; i1<I; i1++){
    dists[i1*I+i1] = 0;
    for (int i2 = i1 +1; i2<I; i2++){
      double d = 0;
      for (int j=0; j<J; j++){
	d += (X[i1*J+j] - X[i2*J+j])*(X[i1*J+j] - X[i2*J+j]);
      }
      dists[i1*I+i2] = d;
      dists[i2*I+i1] = d;
    }
  }
}

void dataset::partitionning(dataset& train, dataset& test, float p){
  int I_test = (int)ceil(I*p);
  int I_train = I - I_test;
  
  test.X = new double[I_test * J];
  test.Y = new int[I_test];
  train.X = new double[I_train * J];
  train.Y = new int[I_train];
  
  vector<int> indexes;
  for (int i=0; i<I; i++){
    indexes.push_back(i);
  }
  random_shuffle(indexes.begin(),indexes.end());
    
  for (int i=0; i<I_test; i++){
    test.Y[i] = Y[indexes[i]];
    for (int j = 0; j < J; j++) {
      test.X[i * J + j] = X[indexes[i] * J + j];
    }
  }
  for (int i=0; i<I_train; i++){
    train.Y[i] = Y[indexes[i+I_test]];
    for (int j = 0; j < J; j++) {
      train.X[i * J + j] = X[indexes[i+I_test] * J + j];
    }
  }
  test.I = I_test;
  test.J = J;
  test.K = K;
  test.name = name;
  
  train.I = I_train;
  train.J = J;
  train.K = K;
  train.compute_mu();
  train.compute_L_hat();
  train.name = name;
}

void dataset::partitionning2(dataset & train, dataset& validation, dataset& test){
  int I_test = (int)ceil(I*0.2);
  int I_validation = I_test;
  int I_train = I - 2*I_test;

  train.X = new double[I_train * J];
  train.Y = new int[I_train];
  
  validation.X = new double[I_validation * J];
  validation.Y = new int[I_validation];
  
  test.X = new double[I_test * J];
  test.Y = new int[I_test];
  
  vector<int> indexes;
  for (int i=0; i<I; i++){
    indexes.push_back(i);
  }
  random_shuffle(indexes.begin(),indexes.end());
    
  for (int i=0; i<I_train; i++){
    train.Y[i] = Y[indexes[i]];
    for (int j = 0; j < J; j++) {
      train.X[i * J + j] = X[indexes[i] * J + j];
    }
  }
  for (int i=0; i<I_validation; i++){
    validation.Y[i] = Y[indexes[i+I_train]];
    for (int j = 0; j < J; j++) {
      validation.X[i * J + j] = X[indexes[i+I_train] * J + j];
    }
  }
  for (int i=0; i<I_test; i++){
    test.Y[i] = Y[indexes[i+I_train+I_validation]];
    for (int j = 0; j < J; j++) {
      test.X[i * J + j] = X[indexes[i+I_train+I_validation] * J + j];
    }
  }
  train.I = I_train;
  train.J = J;
  train.K = K;
  train.compute_mu();
  train.compute_L_hat();
  train.name = name;

  validation.I = I_validation;
  validation.J = J;
  validation.K = K;
  validation.name = name;
  
  test.I = I_test;
  test.J = J;
  test.K = K;
  test.name = name;
}

void dataset::readPartition(int partition_number, dataset& train, dataset& validation, dataset& test){
  fstream file;
  file.open("../TreesAndPartitions/"+name+"/partition"+to_string(partition_number)+".txt", ios::in);
  int Ntrain, Nvalidation, Ntest;
  file >> Ntrain;
  file >> Nvalidation;
  file >> Ntest;
  train = dataset(Ntrain, J, K);
  validation  = dataset(Nvalidation, J, K);
  test = dataset(Ntest, J, K);

  int index;
  for (int i=0; i<Ntrain; i++){
    file >> index;
    for (int j=0; j<J; j++){
      train.X[i*J+j] = X[index*J+j];
    }
    train.Y[i] = Y[index];
  }
  
  for (int i=0; i<Nvalidation; i++){
    file >> index;
    for (int j=0; j<J; j++){
      validation.X[i*J+j] = X[index*J+j];
    }
    validation.Y[i] = Y[index];
  }
  
  for (int i=0; i<Ntest; i++){
    file >> index;
    for (int j=0; j<J; j++){
      test.X[i*J+j] = X[index*J+j];
    }
    test.Y[i] = Y[index];
  }

  train.name = name;
  train.partition = partition_number;
  train.compute_L_hat();
  train.compute_mu();

  validation.name = name;
  validation.partition = partition_number;
  validation.compute_L_hat();
  validation.compute_mu();

  test.name = name;
  test.partition = partition_number;
  test.compute_L_hat();
  test.compute_mu();
}

