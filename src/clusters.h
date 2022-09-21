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
  int nb; // size of the cluster
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
};


struct distCl{
  double dist;
  int id1, id2;

  distCl(){dist = 0.0; id1 = 0; id2=0;}
  distCl(double d, int a, int b) {dist =d; id1 =a; id2 =b;}
};


vector<distCl> initializeDist(dataset& dt);
vector<distCl> initializeWeightedDist(dataset& dt, double* minDist);
/*
void move1(vector<distCl>& ordDist, int index, string direction = "undefined"); // ne finit pas quand c'est trop grand
void move2(vector<distCl>& ordDist, int index); // finit mais c'est le plus lent
void move3(vector<distCl>& ordDist, int index, double newVal); // finit mais ce n'est pas le plus rapide
*/
void move(vector<distCl>& ordDist, int index, double newVal); // finit et c'est le plus rapide !!


class clustering{
 public:
  int J; // dimension
  vector<cluster> clusters;
  map<int,int> clusterOf; // for any data point, we have the cluster it is in
  map<int,int> placeOf; // for any cluster (referenced by its id), we have its position in the vector clusters
  
  clustering(){
    J = 0;
    clusters = {};
    
    map<int,int> clusterOf,
      placeOf;
  }
  clustering(dataset& dt);
  //clustering(clustering& cl);

  clustering clustering_copy();
  
  pair<int,int> mergeClusters(int c1, int c2);
  void breakCluster(int placeOfCluster, vector<vector<int>> div, dataset& initialDt);
  
  void showClusters();
  void showRepartition();
  void showReduction(int I);
  void showPointsTotal();

  void write(string namefile);

  dataset createDt(dataset& initialDt, bool useMedoid=false);
  
};

clustering greedyClustering(dataset& dt, float p);
clustering weightedGreedyClustering(dataset& dt);
clustering homogeneousClustering(dataset& dt, bool useMedoids = false);

