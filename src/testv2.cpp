#include "testv2.h"

void learningWithClustering(string datasetName){
  dataset dt = dataset(datasetName);

  int nbPartitions = 5,
    nbClusterings;

  vector<string> typeOfClustering;

  if (dt.I/2<1000){
    nbClusterings = 7;
    typeOfClustering = {"hierarchicalClustering5", "hierarchicalClustering10","hierarchicalClustering25",
			"betterGreedyClustering", "hierarchicalClusteringGprop",
			"barycenterClustering", "medoidClustering" };
  }
  else{
    nbClusterings = 6;
    typeOfClustering = {"hierarchicalClustering1","hierarchicalClustering5", "hierarchicalClustering10",
			"hierarchicalClustering25","betterGreedyClustering", "hierarchicalClusteringGprop"};
  }


  int DMAX = 4;
  double timeLimitUniv = 1800.0,
    timeLimitMultiv = 300.0;
  
  solution sol;
  time_t now;
  char *date;

  string filename("../results/learningWithClusterings2"+datasetName+".csv");
  ofstream file_out;
  
  file_out.open(filename, ios_base::app);
  file_out << "Dataset,Partition,Model,TimeLimitPerIt,Clustering,%redOfInitCl,D,Time,NbIter,NbIterOpti,";
  file_out << "AvgIterCl,AvgPseudoGap,AvgFinal%Red,errTrain,errTest,bestTree" << endl;
  file_out.close();
 
  for (int p=0; p<nbPartitions; p++){
    dataset dtTrain, dtValidation, dtTest;
    dt.readPartition(p, dtTrain, dtValidation, dtTest);

    training_results_cl tr;

    int  Nmin = (int)floor(dtTrain.I*0.05);

    dtTrain.computeDists();
    
    vector<clustering> cl = {};
    if (dt.I/2>=1000){
      cl.push_back(hierarchicalClustering(dtTrain, 0.01));
    }
    cl.push_back(hierarchicalClustering(dtTrain, 0.05));
    cl.push_back(hierarchicalClustering(dtTrain, 0.1));
    cl.push_back(hierarchicalClustering(dtTrain, 0.25));
    cl.push_back(weightedGreedyClustering(dtTrain));
    cl.push_back(hierarchicalClustering(dtTrain, (double)cl.back().clusters.size()/(double)dtTrain.I));
    if (dt.I/2<1000){
      cl.push_back(homogeneousClustering(dtTrain, false));
      cl.push_back(homogeneousClustering(dtTrain, true));
    }

    for (int s=0; s<2; s++){
      string modelt = "FOCT";
      if (s==1){
	modelt += "H";
      }
      /* // on a deja sans clustering
      now = time(0);
      date = ctime(& now);
      cout << "Partition " << p << " - NoClustering - " << modelt <<" - start : "<< date<<endl;
      
      try {
	training_results trNC= learning(dtTrain, dtValidation, dtTest, baseModel::FOCT, s==0, 4, timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1),Nmin);
	file_out.open(filename, ios_base::app);
	for (int d=0;d<DMAX-1;d++){
	  file_out << datasetName << "," << p << "," << modelt << "," << timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1)<< ",";
	  file_out << "NoClustering" << ",-," << d+2 << "," << tr.time[d] << "," << tr.nbIter[d] << "," << tr.nbOpti[d] << ",-,-,-,";
	  file_out << (double)tr.errTrain[d]/(double)dtTrain.I << "," << (double)tr.errTest[d]/(double)dtTest.I << ",";
	  file_out << tr.bestTree[d] << endl;
	}
	file_out.close();
      } catch(...){
	file_out.open(filename, ios_base::app);
	for (int d=0;d<DMAX-1;d++){
	  file_out << datasetName << "," << p << "," << modelt << "," << timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1)<< ",";
	  file_out << "NoClustering" << ",-," << d+2 << ", ERROR" << endl;
	}
	file_out.close();
	}*/
	
      for (int c=0; c<nbClusterings; c++){
	cl[c].name = typeOfClustering[c];
	
	now = time(0);
	date = ctime(& now);
	cout << "Partition " << p << " - " <<typeOfClustering[c]<< " - " << modelt <<" - start : "<< date<<endl;
	
	try {
	  training_results_cl tr = learningWithClustering(dtTrain, dtValidation, dtTest, cl[c], baseModel::FOCT, s==0, 4,timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1),Nmin);
	  file_out.open(filename, ios_base::app);
	  for (int d=0;d<DMAX-1;d++){
	    file_out << datasetName << "," << p << "," << modelt << "," << timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1)<< ",";
	    file_out << typeOfClustering[c] << "," << (double)cl[c].clusters.size()/(double)dtTrain.I << ",";
	    file_out << d+2 << "," << tr.time[d] << "," << tr.nbIter[d] << "," << tr.nbOpti[d] << "," << tr.nbIterCl[d] << ",";
	    file_out << tr.pseudoGap[d] << "," << tr.finalReduction[d] << "," << (double)tr.errTrain[d]/(double)dtTrain.I << "," << (double)tr.errTest[d]/(double)dtTest.I << ",";
	    file_out << tr.bestTree[d] << endl;
	  }
	  file_out.close();
	} catch(...){
	  file_out.open(filename, ios_base::app);
	  for (int d=0;d<DMAX-1;d++){
	    file_out << datasetName << "," << p << "," << modelt << "," << timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1)<< ",";
	    file_out << typeOfClustering[c] << "," << (double)cl[c].clusters.size()/(double)dtTrain.I << ",";
	    file_out << d+2 << ", ERROR" << endl;
	  }
	  file_out.close();
	}
		
      }
    }
  }
}
