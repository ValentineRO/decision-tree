#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include "utilities.h"
using namespace std;

class dataset{
public :
  int I,J,K;
  string name; // nom du fichier associé au dataset si c'est le cas
  int partition; // numéro de la partition si c'est un dataset de train, test ou validation et -1 si c'est le dataset complet

  double* X; // xij = X[i*J+j]
  int* Y;
  int L_h;
  double mu_min, mu_max;
  double* mu_vect;
    
  double* dists;
  int* repOfLabels;
  int* weights;
  bool weightedPoints;
  int initialI;
    
  dataset();
  dataset(int I_, int J_, int K_);
  dataset(string namefile);
    
  void compute_L_hat();
  void compute_mu();
  void computeDists(); // also compute the repartition of labels

  void partitionning(dataset& train, dataset& test, float p=0.2);
  void partitionning2(dataset& train, dataset& validation, dataset& test);

  void readPartition(int partition_number, dataset& train, dataset& validation, dataset& test);

  void writeDataset(string namefile);
};

