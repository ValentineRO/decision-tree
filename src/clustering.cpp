#include "clustering.h"

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

clustering::clustering(dataset& dt, vector<vector<int>> part){
  J = dt.J;
  clusters = {};
  for (int c=0; c<part.size(); c++){
    int minID = dt.I;
    for (auto &pt: part[c]){
      if (minID > pt){
	minID = pt;
      }
    }
    double* pos = new double[dt.J];
    for (int j=0; j<dt.J; j++){
      pos[j] = dt.X[minID*dt.J+j];
    }
    clusters.push_back(cluster(minID, dt.Y[minID], pos, dt.J));
    for (auto &pt: part[c]){
      if (pt != minID){
	clusters[c].pts.push_back(pt);
	clusterOf[pt] = minID;
      }
    }
    clusters[c].computeBarycenter(dt);
    clusters[c].computeMedoid(dt);
    clusters[c].updateLabel(dt);
    placeOf[minID] = c;
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
  
  for (int j=0; j<J; j++){
    double bar_newJ = 0;
    bar_newJ += clusters[placeOf[c1]].bar[j] * clusters[placeOf[c1]].pts.size();
    bar_newJ += clusters[placeOf[c2]].bar[j] * clusters[placeOf[c2]].pts.size();
    bar_newJ /= clusters[placeOf[c1]].pts.size() + clusters[placeOf[c2]].pts.size();
    clusters[placeOf[c1]].bar[j] = bar_newJ;
  }

  for (int i=0; i<clusters[placeOf[c2]].pts.size(); i++){
    clusters[placeOf[c1]].pts.push_back(clusters[placeOf[c2]].pts[i]);
    clusterOf[clusters[placeOf[c2]].pts[i]] = c1;
  }
  
  clusters[placeOf[c2]] = clusters[clusters.size()-1];
  placeOf[clusters[clusters.size()-1].id] = placeOf[c2];

  clusters.pop_back();
  
  placeOf.erase(c2);
  
  return pair<int,int>(c1,c2);	  
}

void clustering::breakCluster(int placeOfCluster, vector<vector<int>> div, dataset& initialDt){ 
  clusters[placeOfCluster].pts = div[0];
  clusters[placeOfCluster].id = *min_element(div[0].begin(), div[0].end());
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
      }
      clusterOf[div[g][i]] = id;
    }
    clusters[clusters.size()-1].computeBarycenter(initialDt);
    clusters[clusters.size()-1].computeMedoid(initialDt);	       
  }
}

void clustering::showClusters(){
  for (int n=0; n<clusters.size(); n++){
    cout << "Cluster " << n <<endl; 
    for (int i=0; i< clusters[n].pts.size(); i++){
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
    while (rep.size() <= clusters[n].pts.size()){
      rep.push_back(0);
    }

    rep[clusters[n].pts.size()] += 1;
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


void clustering::showPointsTotal(){
  int sum = 0;

  for (int c = 0; c<clusters.size(); c++){
    sum += clusters[c].pts.size();
  }

  cout << "There is a total of " << sum << " points" << endl;
}

float clustering::computeHomogeneity(dataset& dt){
  int nb = 0;
  for (int c=0; c<clusters.size(); c++){
    for (int i=0; i<clusters[c].pts.size(); i++){
      if (dt.Y[clusters[c].pts[i]] == clusters[c].lb){
	nb += 1;
      }
    }
  }
  return (float)nb/dt.I;
}

float clustering::computeExclusion(dataset& dt){
  bool* isIntersected = new bool[clusters.size()];
  for (int c=0; c<clusters.size(); c++){
    isIntersected[c] = false;
  }

  GRBEnv env = GRBEnv();
  int cpt = 0;
  for (int c1=0; c1<clusters.size(); c1++){
    for (int c2=c1+1; c2<clusters.size(); c2++){
      if (not (isIntersected[c1] and isIntersected[c2])){
	GRBModel m = GRBModel(env);

	GRBVar* lambda;
	GRBVar* mu;
	GRBVar* z;

	variable_def setLambda = variable_def("lambda",clusters[c1].pts.size(), false,0,1);
	variable_def setMu = variable_def("mu",clusters[c2].pts.size(), false,0,1);
	variable_def setZ = variable_def("z",dt.J+2, false);

	lambda = m.addVars(setLambda.lb, setLambda.ub, setLambda.coef, setLambda.type, setLambda.names, setLambda.count);
	mu = m.addVars(setMu.lb, setMu.ub, setMu.coef, setMu.type, setMu.names, setMu.count);
	z = m.addVars(setZ.lb, setZ.ub, setZ.coef, setZ.type, setZ.names, setZ.count);

	GRBLinExpr sum_lambda = 0;
	for (int i=0; i<clusters[c1].pts.size(); i++){
	  sum_lambda += lambda[i];
	}
	//m.addConstr(sum_lambda, GRB_EQUAL, 1,"sum_lambda_equals_1");
	m.addConstr(sum_lambda+z[dt.J], GRB_EQUAL, 1,"sum_lambda_equals_1");
	
	GRBLinExpr sum_mu = 0;
	for (int i=0; i<clusters[c2].pts.size(); i++){
	  sum_mu += mu[i];
	}
	m.addConstr(sum_mu+z[dt.J+1], GRB_EQUAL, 1,"sum_mu_equals_1");
	//m.addConstr(sum_mu, GRB_EQUAL, 1,"sum_mu_equals_1");

	for (int j=0; j<dt.J; j++){
	  GRBLinExpr sum_lambdaXj = 0;
	  for (int i=0; i<clusters[c1].pts.size(); i++){
	    sum_lambdaXj += lambda[i]*dt.X[clusters[c1].pts[i]*dt.J+j];
	  }
	  GRBLinExpr sum_muXj = 0;
	  for (int i=0; i<clusters[c2].pts.size(); i++){
	    sum_muXj += mu[i]*dt.X[clusters[c2].pts[i]*dt.J+j];
	  }
	  //m.addConstr(sum_lambdaXj,GRB_EQUAL,sum_muXj,"sum_lambdaXj_equals_sum_muXj_j="+to_string(j));
	  m.addConstr(sum_lambdaXj-sum_muXj+z[j],GRB_EQUAL,0,"sum_lambdaXj_equals_sum_muXj_j="+to_string(j));
	}
	GRBLinExpr obj = 0;
	for (int k=0; k<dt.J+2; k++){
	  obj += z[k];
	}
	m.setObjective(obj,GRB_MINIMIZE);
	//m.write("exclMd"+to_string(cpt)+".lp");
	freopen("intersectingClusters.txt", "a", stdout);
	m.optimize();
	freopen("/dev/tty", "w", stdout);
	remove("intersectingClusters.txt");

	if (m.get(GRB_DoubleAttr_ObjVal) == 0){
	  isIntersected[c1] = true;
	  isIntersected[c2] = true;
	}

	cpt += 1;
      }
    }
  }

  int nb=0;
  for (int c=0; c<clusters.size(); c++){
    if (isIntersected[c]){
      nb += 1;
    }
  }

  return 1-(float)nb/clusters.size();
}

float clustering::computeConsistency(dataset& dt){
  int* closestNeighbour = new int[dt.I];
  for (int i1=0; i1<dt.I; i1++){
    double minDist = -1;
    int cN = -1;
    for (int i2=0; i2<dt.I; i2++){
      if (i2 != i1 and (minDist == -1 or minDist>dt.dists[i1*dt.I+i2])){
	minDist = dt.dists[i1*dt.I+i2];
	cN = i2;
      }
    }
    closestNeighbour[i1] = cN;
  }

  int nb = 0;
  for (int c=0; c<clusters.size(); c++){
    if (clusters[c].pts.size()>1){
      for (auto &pt:clusters[c].pts){
	if (clusters[c].lb != clusters[placeOf[clusterOf[closestNeighbour[pt]]]].lb){
	  nb += 1;
	}
      }
    }
  }

  return 1 - (float)nb/dt.I;
}

float clustering::computeMeanDist(dataset& dt){
  float minDist = 0;
  for (int c=0; c<clusters.size(); c++){
    for (auto &i: clusters[c].pts){
      for (int j=0; j<J; j++){
	minDist += (dt.X[i*J+j]-clusters[c].bar[j])*(dt.X[i*J+j]-clusters[c].bar[j]);
      }
    }
  }
  minDist /= dt.I;
  return minDist;
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


dataset clustering::createDt(dataset& initialDt, bool useMedoid){
  dataset dt = dataset(clusters.size(), J, initialDt.K);

  dt.initialI = initialDt.I;
  dt.weights = new int[dt.I];
  dt.weightedPoints = true;

  for (int i=0; i<dt.I; i++){
    dt.Y[i] = clusters[i].lb;
    dt.weights[i] = clusters[i].pts.size();
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

clustering hierarchicalClustering(dataset& dt, float p, bool useLabel){
  clustering C = clustering(dt);
  vector<distCl> ordDist = initializeDist(dt); // here d.id1 and d.id2 are the data points id

  while (((float)C.clusters.size()/(float)dt.I > p) and (ordDist.size() >0)){
    bool canMerge = (C.clusterOf[ordDist.back().id1] != C.clusterOf[ordDist.back().id2]);
    if (useLabel){
      canMerge = canMerge and dt.Y[ordDist.back().id1] == dt.Y[ordDist.back().id2];
    }    
    if (canMerge){
      int c1 = C.clusterOf[ordDist.back().id1],
	c2 = C.clusterOf[ordDist.back().id2];
      pair<int,int> c = C.mergeClusters(c1,c2);
      C.clusters[C.placeOf[c.first]].computeMedoid(dt);
      if (not useLabel){
	C.clusters[C.placeOf[c.first]].updateLabel(dt);
      }
    }
    ordDist.pop_back();
  }
  return C;
}

clustering auxKMeans(dataset& dt, float p, bool useLabel){
  random_device rd;
  default_random_engine eng(rd());
  uniform_real_distribution<float> distr(0.0, 1.0);
  
  int nbCl = (int)ceil(dt.I*p);
  float* posCenters = new float[nbCl*dt.J];
  for (int c=0; c<nbCl; c++){
    for (int j=0; j<dt.J; j++){
      posCenters[c*dt.J+j] = distr(eng);
    }
  }

  vector<int> labels;
  if (useLabel){
    for (int k=0; k<dt.K; k++){
      int nbClForK = round((float) dt.repOfLabels[k] / dt.I * nbCl);
      for (int it=0; it < nbClForK; nbClForK){
	labels.push_back(k);
      }
    }
    random_shuffle(labels.begin(), labels.end());
  }

  int* associatedCenter = new int[dt.I];
  for (int i=0; i<dt.I; i++){
    associatedCenter[i] = -1;
  }
  vector<vector<int>> tempClusters = {};
  
  bool asChanged = true;

  while (asChanged){
    asChanged = false;
    if (associatedCenter[0] != -1){// we dont want to be computing new center for the first step of the algorithm
      for (int c=0; c<nbCl; c++){
	for (int j=0; j<dt.J; j++){
	  posCenters[c*dt.J+j] = 0;
	  for (auto &i: tempClusters[c]){
	    posCenters[c*dt.J+j] += dt.X[i*dt.J+j];
	  }
	  posCenters[c*dt.J+j] /= tempClusters[c].size();
	}
      }
    }
    tempClusters = {};
    for (int c=0; c<nbCl; c++){
      tempClusters.push_back({});
    }
    for (int i=0; i<dt.I; i++){
      float minDist = -1;
      int newCenter = -1;
      for (int c=0; c<nbCl; c++){
	bool consideringClusterC = true;
	if (useLabel){
	  consideringClusterC = (labels[c] == dt.Y[i]);
	}
	if (consideringClusterC){
	  float dist = 0;
	  for (int j=0; j<dt.J; j++){
	    dist += (dt.X[i*dt.J+j] - posCenters[c*dt.J+j])*(dt.X[i*dt.J+j] - posCenters[c*dt.J+j]);
	  }
	  if (minDist == -1 or minDist > dist){
	    minDist = dist;
	    newCenter = c;
	  }
	}
      }
      if (newCenter != associatedCenter[i]){
	asChanged = true;
	associatedCenter[i] = newCenter;
      }
      tempClusters[associatedCenter[i]].push_back(i);
    }
  }

  return clustering(dt, tempClusters);
}

clustering kMeansClustering(dataset& dt, float p, bool useLabel){
  clustering bestKMeans;
  float score = 0;
  for (int cl=0; cl<10; cl++){
    clustering newKMeans = auxKMeans(dt,p,useLabel);
    
  }

  return bestKMeans;
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
        C.clusters[C.placeOf[c.first]].computeMedoid(dt);
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
      
      C.clusters[C.placeOf[c.first]].computeMedoid(dt);
      
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

clustering optimalClustering(dataset& dt, int maxCl, float H){
  int CMAX = maxCl;
  if (maxCl == -1){
    CMAX = dt.I;
  }
  
  GRBEnv env = GRBEnv();
  GRBModel md = GRBModel(env);

  int* cN = new int[dt.I]; // closest neighbor of i
  bool* B = new bool[dt.I]; // is the closest neighbor of i of identical label
  float maxDist = 0;

  for (int i=0; i<dt.I; i++){
    float minDist = -1;
    for (int j=0; j<dt.I; j++){
      if ((minDist == -1 or minDist>dt.dists[i*dt.I+j]) and i != j){
	cN[i] = j;
	minDist = dt.dists[i*dt.I+j];
      }
      if (dt.dists[i*dt.I+j] > maxDist){
	maxDist = dt.dists[i*dt.I+j];
      }
    }
    B[i] = dt.Y[i] == dt.Y[cN[i]];
  }

  GRBVar* z = md.addVars(dt.I * CMAX, GRB_BINARY); // z_ic = z[i*CMAX+c] : i belongs to c
  GRBVar* r = md.addVars(dt.I * CMAX, GRB_BINARY); // r_ic = r[i*CMAX+c] : i represents c
  GRBVar* D = md.addVars(dt.I, GRB_CONTINUOUS); // distance between i and its rep c
  GRBVar* DMAX = md.addVars(CMAX, GRB_CONTINUOUS); // max distance between the rep of c and an element of c
  GRBVar* Nc = md.addVars(CMAX, GRB_CONTINUOUS); // number of element in cluster c
  GRBVar* Nkc = md.addVars(CMAX * dt.K, GRB_CONTINUOUS); // N_kc = Nkc[c*dt.K+k] : number of element of label k in cluster c
  GRBVar* cc = md.addVars(CMAX * dt.K, GRB_BINARY); // c_kc = c[c*dt.K+k] : c is labelled k
  GRBVar* ci;

  if (H < 1){
    ci = md.addVars(dt.I * dt.K, GRB_BINARY); // c_ki = ci[i*dt.K+k] : i is labelled k
  }

  string constraint_name;
  for (int i=0; i<dt.I; i++){ // sum_c zic <= 1
    GRBLinExpr sum_zic = 0;
    for (int c=0; c<CMAX; c++){
      sum_zic += z[i*CMAX+c];
    }
    constraint_name = "sum_zic_leq_1_i="+to_string(i);
    int secondPart = 0;
    if (H<1 or B[i]){
      secondPart = 1;
    }
    md.addConstr(sum_zic,GRB_LESS_EQUAL,secondPart,constraint_name);
  }
  for (int c=0; c<CMAX; c++){ // sum_i ric <= 1
    GRBLinExpr sum_ric = 0;
    for (int i=0; i<dt.I; i++){
      sum_ric += r[i*CMAX+c];
    }
    constraint_name = "sum_ric_leq_1_c="+to_string(c);
    md.addConstr(sum_ric,GRB_LESS_EQUAL,1,constraint_name);
  }
  for (int i=0; i<dt.I; i++){
    for (int c=0; c<CMAX; c++){
      // r_ic <= z_ic
      constraint_name = "ric_leq_zic_i="+to_string(i)+"_c="+to_string(c);
      md.addConstr(r[i*CMAX+c],GRB_LESS_EQUAL,z[i*CMAX+c],constraint_name);
      
      // Di >= sum_i' r_i'c*dist[i,i'] - maxDist*zic
      GRBLinExpr sum_ripcXdist = 0;
      for (int ip=0; ip<dt.I; ip++){
	sum_ripcXdist += dt.dists[i*dt.I+ip]*r[ip*CMAX+c];
      }
      constraint_name = "Di_geq_sum_ip_ripc*dist[ip,i]_-_Mzic_i="+to_string(i)+"_c="+to_string(c);
      md.addConstr(D[i], GRB_GREATER_EQUAL,sum_ripcXdist-maxDist*z[i*CMAX+c],constraint_name);
      
      // DMAXc >= Di - maxDist*zic
      constraint_name = "DMAXc_geq_Di_-_Mzic_i="+to_string(i)+"_c="+to_string(c);
      md.addConstr(DMAX[c], GRB_GREATER_EQUAL,D[i]-maxDist*z[i*CMAX+c],constraint_name);

      for (int ip=0; ip<dt.I; ip++){
	if (ip!=i){// DMAXc <= dist[i,ip] + maxDist(1-r_ic+z_ipc)
	  constraint_name = "DMAXc_leq_dist[i,ip]_+_maxDist*(1_-_ric+zipc)_i="+to_string(i)+"_c="+to_string(c)+"_ip="+to_string(ip);
	  md.addConstr(DMAX[c], GRB_LESS_EQUAL, dt.dists[i*dt.I+ip] + maxDist*(1-r[i*CMAX+c]+z[i*CMAX+c]), constraint_name);
	}
      }
    }
  }
  for (int c=0; c<CMAX; c++){
    // Nc = sum_zic
    GRBLinExpr sum_zic = 0;
    for (int i=0; i<dt.I; i++){
      sum_zic += z[i*CMAX+c];
    }
    constraint_name = "Nc=sum_zic_c="+to_string(c);
    md.addConstr(Nc[c],GRB_EQUAL, sum_zic,constraint_name);
    
    for (int k=0; k<dt.K; k++){
      // Nkc = sum_zic*1[Y[i]=k]
      GRBLinExpr sum_zic_K = 0;
      for (int i=0; i<dt.I; i++){
	if (dt.Y[i] == k){
	  sum_zic_K += z[i*CMAX+c];
	}
      }
      constraint_name = "Nkc=sum_zic*1[Y[i]=k]_c="+to_string(c)+"_k="+to_string(k);
      md.addConstr(Nc[c],GRB_EQUAL, sum_zic,constraint_name);

      // Nkc >= H*Nc - |I|(1-c_kc)
      constraint_name = "Nkc_geq_H*Nc-I*(1-c_kc)_c="+to_string(c)+"_k="+to_string(k);
      md.addConstr(Nkc[c*dt.K+k],GRB_GREATER_EQUAL, H*Nc[c] - dt.I*(1-cc[c*dt.K+k]),constraint_name);
    }
    // sum c_kc = 1
    GRBLinExpr sum_ckc = 0;
    for (int k=0; k<dt.K; k++){
      sum_ckc += cc[c*dt.K+k];
    }
    constraint_name = "sum_ckc_eq_1_c="+to_string(c);
    md.addConstr(sum_ckc,GRB_EQUAL,1,constraint_name);
  }
  if (H<1){
    for (int i=0; i<dt.I; i++){
      // ci_ki = ci_kCN[i]
      for (int k=0; k<dt.K; k++){ 
	constraint_name = "ci_k,i_=_ci_k,cN[i]_i="+to_string(i)+"_k="+to_string(k);
	md.addConstr(ci[i*dt.K+k], GRB_EQUAL,ci[cN[i]*dt.K+k],constraint_name);
      }
      // sum_ci_ki = 1
      GRBLinExpr sum_ci_ki = 0;
      for (int k=0; k<dt.K; k++){
	sum_ci_ki += ci[i*dt.K+k];
      }
      constraint_name = "sum_ci_ki_=_1_i="+to_string(i);
      md.addConstr(sum_ci_ki,GRB_EQUAL,1,constraint_name);
      // ci_ki >= c_kc - 2*(1-zic)
      for (int c=0; c<CMAX; c++){
	for (int k=0; k<dt.K; k++){
	  constraint_name = "ci_ki_geq_c_kc_-_2*(1-z_ic)_i="+to_string(i)+"_c="+to_string(c)+"_k="+to_string(k);
	  md.addConstr(ci[i*dt.K+k],GRB_GREATER_EQUAL,cc[c*dt.K+k] - 2*(1-z[i*CMAX+c]),constraint_name);
	}
      }
    }
  }

  GRBLinExpr obj = 0; // sum_i(1-sum_c zic) + sum_i,c ric
  for (int i=0; i<dt.I; i++){
    obj += 1;
    for (int c=0; c<CMAX; c++){
      obj -= z[i*CMAX+c];
      obj += r[i*CMAX+c];
    }
  }
  md.setObjective(obj, GRB_MINIMIZE);
  md.optimize();

  vector<vector<int>> partition;
  int cpt = 0;
  for (int c=0; c<CMAX; c++){
    partition.push_back({});
    for (int i=0; i<dt.I; i++){
      if (z[i*CMAX+c].get(GRB_DoubleAttr_X) > 0.9){
	partition[cpt].push_back(i);
      }
    }
    if (partition[cpt].size() == 0){
      partition.pop_back();
    }
    else{
      cpt += 1;
    }
  }
  for (int i=0; i<dt.I; i++){
    bool isNotInCl = true;
    for (int c=0; c<CMAX; c++){
      if (z[i*CMAX+c].get(GRB_DoubleAttr_X) > 0.9){
	isNotInCl = false;
	break;
      }
    }
    if (isNotInCl){
      partition.push_back({i});
    }
  }

  return clustering(dt, partition);
}
