#include "clusters.h"

cluster::cluster(int id1, int label, double* pos1, int dim_pos){
  id = id1;
  m = id1;
  nb = 1;
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
  for (int i=0; i<nb; i++){
    double d = 0;
    for (int j=0; j<J; j++){
      d += (bar[j] - dt.X[pts[i]*J+j])*(bar[j] - dt.X[pts[i]*J+j]);
    }

    if (d<distMin){
      distMin = d;
      ind = pts[i];
    }
  }

  m = ind;
}

cluster cluster::cluster_copy(){
  cluster cl = cluster(id, lb, bar, J);

  for (int i=1; i<pts.size(); i++){
    cl.pts.push_back(pts[i]);
  }
  cl.nb = nb;
  cl.m = m;

  return cl;
}

vector<distCl> initializeDist(dataset& dt){
  vector<distCl> v;

  for (int i1=0; i1<dt.I; i1++){
    for (int i2= i1+1; i2<dt.I; i2++){
      distCl dC;
      dC.dist = dt.dists[i1*dt.I+i2];
      dC.id1 = i1;
      dC.id2 = i2;
      v.push_back(dC);
    }
  }

  sort(v.begin(), v.end(), [](const distCl& d1, const distCl& d2) { return d1.dist > d2.dist; });

  return v;
}

vector<distCl> initializeWeightedDist(dataset& dt, double* minDist){
  for (int i1=0; i1<dt.I; i1++){
    minDist[i1] = 2*dt.J;
    for (int i2 = 0; i2<dt.I; i2++){
      if (dt.Y[i1] != dt.Y[i2] and dt.dists[i1*dt.I+i2] < minDist[i1]){
	minDist[i1] = dt.dists[i1*dt.I+i2];
      }
    }
  }

  vector<distCl> v;

  for (int i1=0; i1<dt.I; i1++){
    for (int i2= i1+1; i2<dt.I; i2++){
      distCl dC;
      dC.dist = dt.dists[i1*dt.I+i2]/min(minDist[i1],minDist[i2]);
      dC.id1 = i1;
      dC.id2 = i2;

      if (dC.dist<1){
	v.push_back(dC);
      }
    }
  }

  sort(v.begin(), v.end(), [](const distCl& d1, const distCl& d2) { return d1.dist > d2.dist; });

  return v;
}

/*
void move1(vector<distCl>& ordDist, int index, string direction){ // recursive function
  if (direction == "-"){
    if (index > 0){ // if we were going left and the index is 0, we cant go further
      if (ordDist[index-1].dist< ordDist[index].dist){ // in that case we want to keep going
	// we swap index-1 and index
	distCl temp = ordDist[index];
	ordDist[index] = ordDist[index-1];
	ordDist[index-1] = temp;
	// and we go to the next step
	move(ordDist, index-1, "-");
      }
    }
  }
  else if (direction == "+"){
    if (index < ordDist.size()-1){
      if (ordDist[index+1].dist> ordDist[index].dist){ // in that case we want to keep going
	// we swap index+1 and index
	distCl temp = ordDist[index];
	ordDist[index] = ordDist[index+1];
	ordDist[index+1] = temp;
	// and we go to the next step
	move(ordDist, index+1, "+");
      }
    }
  }
  else{
    if (index == 0){
      move(ordDist, index, "+");
    }
    else if (index == ordDist.size()-1){
      move(ordDist, index, "-");
    }
    else{
      if (ordDist[index-1].dist< ordDist[index].dist){ // in that case we want to keep going
	// we swap index-1 and index
	distCl temp = ordDist[index];
	ordDist[index] = ordDist[index-1];
	ordDist[index-1] = temp;
	// and we go to the next step
	move(ordDist, index-1, "-");
      }
      if (ordDist[index+1].dist> ordDist[index].dist){ // in that case we want to keep going
	// we swap index+1 and index
	distCl temp = ordDist[index];
	ordDist[index] = ordDist[index+1];
	ordDist[index+1] = temp;
	// and we go to the next step
	move(ordDist, index+1, "+");
      }
    }
  }
}

void move2(vector<distCl>& ordDist, int index){
  string direction;
  if (index == 0){
    direction = "+";
  }
  else if (index == ordDist.size()-1){
    direction = "-";
  }
  else{
    if (ordDist[index-1].dist< ordDist[index].dist){
      direction = "-";
    }
    else if (ordDist[index+1].dist> ordDist[index].dist){
      direction = "+";
    }
    else{
      direction = "none";
    }
  }

  int i= index;
  
  if (direction == "+"){
    while(i<ordDist.size()-1 and ordDist[index+1].dist> ordDist[index].dist){
      distCl temp = ordDist[i];
      ordDist[i] = ordDist[i+1];
      ordDist[i+1] = temp;
      i += 1;
    }
  }

  if (direction == "-"){
    while (i>0 and ordDist[index-1].dist< ordDist[index].dist){
      distCl temp = ordDist[i];
      ordDist[i] = ordDist[i-1];
      ordDist[i-1] = temp;
      i -= 1;
    }
  }
}

void move3(vector<distCl>& ordDist, int index, double newVal){
  int pos = -1;
  if (newVal >= ordDist[0].dist){
    pos = 0;
  }
  else if (newVal < ordDist.back().dist){
    pos = ordDist.size();
  }
  else{
    int a=0,
      b= ordDist.size()-1;
    int c = (a+b)/2;

    while (newVal > ordDist[c].dist or newVal < ordDist[c+1].dist) {
      if (newVal > ordDist[c].dist){
	b = c;
	c= (a+b)/2;
	if (ordDist[b].dist == newVal){
	  pos = b;
	  break;
	}
      }
      else{
	a = c;
	c= (a+b)/2;
	if (ordDist[a].dist == newVal){
	  pos = a;
	  break;
	}
      }
      
      if (b-a == 1){
	break;
      }
    }

    if (pos == -1){
      if (ordDist[c].dist == newVal){
	pos = c;
      }
      else{
	pos = c+1;
      }
    }
    
  }

  if (pos == index){
    ordDist[index].dist = newVal;
  }
  else if (pos<index){
    distCl temp = ordDist[index];
    temp.dist = newVal;
    ordDist.erase(ordDist.begin()+index);
    ordDist.insert(ordDist.begin()+pos,temp);
  }
  else{
    distCl temp = ordDist[index];
    temp.dist = newVal;
    pos -= 1;
    ordDist.erase(ordDist.begin()+index);
    ordDist.insert(ordDist.begin()+pos,temp);
  }
}

*/

void move(vector<distCl>& ordDist, int index, double newVal){
  int pos = -1;
  if (newVal >= ordDist[0].dist){
    pos = 0;
  }
  else if (newVal < ordDist.back().dist){
    pos = ordDist.size();
  }
  else{
    int a=0,
      b= ordDist.size()-1;
    int c = (a+b)/2;

    while (newVal > ordDist[c].dist or newVal < ordDist[c+1].dist) {
      if (newVal > ordDist[c].dist){
	b = c;
	c= (a+b)/2;
	if (ordDist[b].dist == newVal){
	  pos = b;
	  break;
	}
      }
      else{
	a = c;
	c= (a+b)/2;
	if (ordDist[a].dist == newVal){
	  pos = a;
	  break;
	}
      }
      
      if (b-a == 1){
	break;
      }
    }

    if (pos == -1){
      if (ordDist[c].dist == newVal){
	pos = c;
      }
      else{
	pos = c+1;
      }
    }
    
  }

  if (pos == index){
    ordDist[index].dist = newVal;
  }
  else if (pos<index){
    distCl temp = ordDist[index];
    temp.dist = newVal;
    for (int k=index; k>pos; k--){
      ordDist[k] = ordDist[k-1];
    }
    ordDist[pos] = temp;
  }
  else{
    distCl temp = ordDist[index];
    temp.dist = newVal;
    for (int k=index; k<pos-1;k++){
      ordDist[k] = ordDist[k+1];
    }
    ordDist[pos-1] = temp;
  }
}

clustering::clustering(dataset& dt){
  J = dt.J;
  clusters = {};

  for (int i=0; i<dt.I; i++){
    double* pos = new double[dt.J];
    for (int j=0; j<dt.J; j++){
      pos[j] = dt.X[i*dt.J+j];
    }
    clusters.push_back(cluster(i, dt.Y[i], pos, dt.J));
    clusterOf[i]=i;
    placeOf[i]=i;
  }
}


clustering clustering::clustering_copy(){
  clustering cl = clustering();
  
  cl.J = J;

  cl.clusterOf = clusterOf;
  cl.placeOf = placeOf;

  cl.clusters = {};
  for (int i=0; i<clusters.size();i++){
    cl.clusters.push_back(clusters[i].cluster_copy());
  }

  return cl;
}

pair<int,int> clustering::mergeClusters(int c1, int c2){ // note that c1 and c2 are ids of clusters
  if (c1 > c2){ // we want to keep the smallest id
    int temp =  c1;
    c1 = c2;
    c2 = temp;
  }  
  for (int i=0; i<clusters[placeOf[c2]].nb; i++){
    clusters[placeOf[c1]].pts.push_back(clusters[placeOf[c2]].pts[i]);
    clusterOf[clusters[placeOf[c2]].pts[i]] = c1;
  }

  for (int j=0; j<J; j++){
    clusters[placeOf[c1]].bar[j] = (clusters[placeOf[c1]].bar[j]*clusters[placeOf[c1]].nb + clusters[placeOf[c2]].bar[j]*clusters[placeOf[c2]].nb)/(clusters[placeOf[c1]].nb+clusters[placeOf[c2]].nb);
  }

  clusters[placeOf[c1]].nb += clusters[placeOf[c2]].nb;

  clusters[placeOf[c2]] = clusters[clusters.size()-1];
  placeOf[clusters[clusters.size()-1].id] = placeOf[c2];

  clusters.pop_back();
  
  placeOf.erase(c2);
  
  return pair<int,int>(c1,c2);	  
}

void clustering::breakCluster(int placeOfCluster, vector<vector<int>> div, dataset& initialDt){ 
  clusters[placeOfCluster].pts = div[0];
  clusters[placeOfCluster].id = *min_element(div[0].begin(), div[0].end());
  clusters[placeOfCluster].nb = div[0].size();
  clusters[placeOfCluster].computeBarycenter(initialDt);
  clusters[placeOfCluster].computeMedoid(initialDt);

  for (int g=1; g<div.size(); g++){
    int id = *min_element(div[g].begin(), div[g].end());
    double* pos = new double[initialDt.J];
    clusters.push_back(cluster(id, initialDt.Y[id], pos, initialDt.J));
    placeOf[id] = clusters.size()-1;
    for (int i=0; i<div[g].size(); i++){
      if (div[g][i] != id){
	clusters[clusters.size()-1].pts.push_back(div[g][i]);
	clusters[clusters.size()-1].nb += 1;
      }
      clusterOf[div[g][i]] = id;
    }
    clusters[clusters.size()-1].computeBarycenter(initialDt);
    clusters[clusters.size()-1].computeMedoid(initialDt);	       
  }
}

void clustering::showClusters(){
  for (int n=0; n<clusters.size(); n++){
    cout << "Cluster " << n << endl; 
    for (int i=0; i< clusters[n].nb; i++){
      cout << clusters[n].pts[i] << "  ";
    }
    cout << endl;
  }
  cout << endl;
}

void clustering::showRepartition(){
  vector<int> rep;
  rep.push_back(0);
  
  for (int n=0; n<clusters.size(); n++){
    while (rep.size() <= clusters[n].nb){
      rep.push_back(0);
    }

    rep[clusters[n].nb] += 1;
  }
  
  for (int i=1; i<rep.size(); i++){
    if (rep[i] != 0){
      cout << "Clusters of size " << i << ": " << rep[i] << endl;
    }
  }

  cout << endl;
}

void clustering::showReduction(int I){
  cout << 100*(float)clusters.size()/(float)I << "% of reduction" << endl;
}

void clustering::write(string namefile){
  fstream file;
  file.open(namefile,ios::out);

  for (int i=0; i<clusters.size(); i++){
    file << clusters[i].pts.size() << "\n";
    for (int j=0; j<clusters[i].pts.size(); j++){
      file << clusters[i].pts[j] << "\n";
    }
  }
  file << "END";
  file.close();
}

void clustering::showPointsTotal(){
  int sum = 0;

  for (int c = 0; c<clusters.size(); c++){
    sum += clusters[c].pts.size();
  }

  cout << "There is a total of " << sum << " points" << endl;
}

dataset clustering::createDt(dataset& initialDt, bool useMedoid){
  dataset dt = dataset(clusters.size(), J, initialDt.K);

  dt.initialI = initialDt.I;
  dt.weights = new int[dt.I];

  for (int i=0; i<dt.I; i++){
    dt.Y[i] = clusters[i].lb;
    dt.weights[i] = clusters[i].nb;
    for (int j=0; j<J; j++){
      if (useMedoid){
	dt.X[i*J+j]= initialDt.X[clusters[i].m*J+j];
      }
      else{
	dt.X[i*J+j] = clusters[i].bar[j];
      }
    }
  }

  dt.compute_mu();
  dt.compute_L_hat();

  return dt;
}

clustering greedyClustering(dataset& dt, float p){
  clustering C = clustering(dt);
  vector<distCl> ordDist = initializeDist(dt); // here d.id1 and d.id2 are the data points id

  while (((float)C.clusters.size()/(float)dt.I > p) and (ordDist.size() >0)){
    if ((dt.Y[ordDist.back().id1] == dt.Y[ordDist.back().id2]) and (C.clusterOf[ordDist.back().id1] != C.clusterOf[ordDist.back().id2])){
      int c1 = C.clusterOf[ordDist.back().id1],
	c2 = C.clusterOf[ordDist.back().id2];
      pair<int,int> c = C.mergeClusters(c1,c2);
    }
    ordDist.pop_back();
  }

  return C;
}

clustering weightedGreedyClustering(dataset& dt){
  clustering C = clustering(dt);

  double* minDist = new double[dt.I]; // identify the closest neighbor of different class

  vector<distCl> ordDist = initializeWeightedDist(dt, minDist); // here d.id1 and d.id2 are the data points id

  while (ordDist.size()>0){
    int i1 = ordDist.back().id1,
      i2 = ordDist.back().id2;
    int indexMedC1 = C.clusters[C.placeOf[C.clusterOf[i1]]].m,
      indexMedC2 = C.clusters[C.placeOf[C.clusterOf[i2]]].m;
    if (minDist[indexMedC1] != 0 and minDist[indexMedC2] != 0){
      if (dt.dists[indexMedC1*dt.I+indexMedC2]/minDist[indexMedC1] < 1 and dt.dists[indexMedC1*dt.I+indexMedC2]/minDist[indexMedC2] < 1 and indexMedC1 != indexMedC2){
	int c1 = C.clusterOf[i1],
	  c2 = C.clusterOf[i2];
	pair<int,int> c = C.mergeClusters(c1,c2);
	C.clusters[c.first].computeMedoid(dt);
      }
    }
    ordDist.pop_back();
  }
  

  return C;
}

clustering homogeneousClustering(dataset& dt, bool useMedoids){
  clustering C = clustering(dt);
  vector<distCl> ordDist = initializeDist(dt); // here d.id1 and d.id2 are the id of the cluster
  
  int mergeableClusters = dt.I;
  bool* isMergeable = new bool[dt.I];

  for (int i=0; i<dt.I; i++){
    isMergeable[i] = true;
  }
  
  while(mergeableClusters > 1){
    if (dt.Y[ordDist.back().id1] != dt.Y[ordDist.back().id2]){
      int deletedClusters = (int)isMergeable[ordDist.back().id1] + (int)isMergeable[ordDist.back().id2];
      mergeableClusters -= deletedClusters;
      
      isMergeable[ordDist.back().id1] = false;
      isMergeable[ordDist.back().id2] = false;
      ordDist.pop_back();
    }
    else if (isMergeable[ordDist.back().id1] and isMergeable[ordDist.back().id2]){
      pair<int,int> c = C.mergeClusters(ordDist.back().id1,ordDist.back().id2);
      ordDist.pop_back();
      
      if (useMedoids){
	C.clusters[c.first].computeMedoid(dt);
      }

      // let's remove distances linked to c2 ie c.second
      int cpt = ordDist.size()-1;

      while (cpt >= 0){
	if ((ordDist[cpt].id1 == c.second) or (ordDist[cpt].id2 == c.second)){
	  ordDist.erase(ordDist.begin()+cpt);
	}
	cpt -= 1;
      }

      isMergeable[c.second] = false;
      mergeableClusters -= 1;

      // let's change distances linked to c1 ie c.first
      for (int i=0; i<ordDist.size(); i++){
	if ((ordDist[i].id1 == c.first) or (ordDist[i].id2 == c.first)){
	  double d=0;
	  if (useMedoids){
	    for (int j=0; j<dt.J; j++){
	      d += pow(dt.X[C.clusters[C.placeOf[ordDist[i].id1]].m*C.J+j] - dt.X[C.clusters[C.placeOf[ordDist[i].id2]].m*C.J+j],2);
	    }
	  }
	  else{// otherwise we use barycenters
	    for (int j=0; j<dt.J; j++){
	      d += pow(C.clusters[C.placeOf[ordDist[i].id1]].bar[j] - C.clusters[C.placeOf[ordDist[i].id2]].bar[j],2);
	    }
	  }
	  move(ordDist,i,d);
	}
      }
    }
    else{
      ordDist.pop_back();
    }
  }
  return C;
}
