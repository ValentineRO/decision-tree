#pragma once
#include <random>
#include "clusters.h"

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
  string name;
  int J; // dimension
  vector<cluster> clusters;
  map<int,int> clusterOf; // for any data point, we have the cluster it is in
  map<int,int> placeOf; // for any cluster (referenced by its id), we have its position in the vector clusters;
  dataset* refDt;
  
  clustering(){
    J = 0;
    clusters = {};
    
    map<int,int> clusterOf,
      placeOf;
  }
  clustering(dataset& dt);
  clustering(dataset& dt, vector<vector<int>> part);
  clustering(dataset& dt, string namefile, string shortName);

  clustering clustering_copy();
  
  pair<int,int> mergeClusters(int c1, int c2);
  void breakCluster(int placeOfCluster, vector<vector<int>> div, dataset& initialDt, bool variableRep=false);
  
  void showClusters();
  void showRepartition();
  void showReduction(int I);
  void showPointsTotal();

  float computeHomogeneity(dataset& dt);
  float computeExclusion(dataset& dt);
  float computeConsistency(dataset& dt);
  float computeMeanDist(dataset& dt);

  void write(string namefile);

  dataset createDt(dataset& initialDt, bool useMedoid=false, bool useRep=false);
};

clustering hierarchicalClustering(dataset& dt, float p, bool useLabel=true);
clustering pythonHierarchicalClustering(dataset& dt, float p);
clustering auxKMeans(dataset& dt, float p, bool useLabel=true);
clustering kMeansClustering(dataset& dt, float p, bool useLabel=true);
clustering weightedGreedyClustering(dataset& dt);
clustering homogeneousClustering(dataset& dt, bool useMedoids = false);
clustering optimalClustering(dataset& dt, int maxCl=-1, float H=1);
