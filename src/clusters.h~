#pragma once
#include "model_utilities.h"
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
using namespace std;

class cluster{
 public:
  int id; // id of the cluster, it will be the smallest index
  vector<int> pts; // list of the indexes of the data point in the cluster
  int m; // index of the medoid
  double* bar; // position of the barycenter
  int J;
  int lb;

  cluster();
  cluster(int id1, int label, double* pos1, int dim_pos);
  //cluster(cluster& c);

  cluster cluster_copy();
  
  void computeBarycenter(dataset& dt);
  void computeMedoid(dataset& dt);
  void updateLabel(dataset& dt);
};


struct distCl{
  double dist;
  int id1, id2;

  distCl(){dist = 0.0; id1 = 0; id2=0;}
  distCl(double d, int a, int b) {dist =d; id1 =a; id2 =b;}
};
