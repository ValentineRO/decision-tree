#include "clusters.h"

cluster::cluster(int id1, int label, double* pos1, int dim_pos){
  id = id1;
  m = id1;
  J = dim_pos;

  lb = label;

  pts = {id1};
  bar = new double[dim_pos];
  for (int j=0; j<dim_pos; j++){
    bar[j] = pos1[j];
  }
}

void cluster::computeBarycenter(dataset& dt){
  for (int j=0; j<dt.J; j++){
    bar[j] = 0;
    for (int i=0; i<pts.size(); i++){
      bar[j] += dt.X[pts[i]*dt.J+j];
    }
    bar[j] /= pts.size();
  }
}

void cluster::computeMedoid(dataset& dt){
  double distMin = 10;
  int ind = -1;
  for (auto i: pts){
    double d = 0;
    for (int j=0; j<J; j++){
      d += (bar[j] - dt.X[i*J+j])*(bar[j] - dt.X[i*J+j]);
    }

    if (d<distMin){
      distMin = d;
      ind = i;
    }
  }

  m = ind;
}

void cluster::updateLabel(dataset& dt){
  int* nbPerLabel = new int[dt.K];
  for (int k=0; k<dt.K; k++){
    nbPerLabel[k] = 0;
  }

  for (auto &pt: pts){
    nbPerLabel[dt.Y[pt]] += 1;
  }

  int maxNb = 0;
  for (int k=0; k<dt.K; k++){
    if (nbPerLabel[k]>maxNb){
      lb = k;
      maxNb = nbPerLabel[k];
    }
  }
}

cluster cluster::cluster_copy(){
  cluster cl = cluster(id, lb, bar, J);

  for (auto i: pts){
    if (i != id){
      cl.pts.push_back(i);
    }
  }
  cl.m = m;

  return cl;
}
