#include "testv2.h"

void learningWithClustering(string datasetName){
  dataset dt = dataset(datasetName);

  int nbClusterings = 8,
    nbPartitions = 5;

  string typeOfClustering[8] = {"hierarchicalClustering5", "hierarchicalClustering10","hierarchicalClustering25",
				"hierarchicalClustering50", "barycenterClustering", "medoidClustering", "betterGreedyClustering",
				"hierarchicalClusteringGprop"};

  int DMAX = 4;
  double timeLimitUniv = 1800.0,
    timeLimitMultiv = 900.0;
  
  solution sol;
  time_t now;
  char *date;

  string filename("../results/learningWithClusterings"+datasetName+".csv");
  ofstream file_out;
  
  file_out.open(filename, ios_base::app);
  file_out << "Dataset,Partition,Model,TimeLimitPerIt,Clustering,SelectingStrat,%redOfInitCl,D,Time,NbIter,NbIterOpti,";
  file_out << "AvgIterCl,AvgPseudoGap,AvgFinal%Red,errTrain,errTest,bestTree" << endl;
  file_out.close();
 
  for (int p=0; p<nbPartitions; p++){
    dataset dtTrain, dtValidation, dtTest;
    dt.readPartition(p, dtTrain, dtValidation, dtTest);

    training_results_cl tr;

    int  Nmin = (int)floor(dtTrain.I*0.05);

    dtTrain.computeDists();
    
    clustering* cl = new clustering[8];
    cl[0] = hierarchicalClustering(dtTrain, 0.05);
    cl[1] = hierarchicalClustering(dtTrain, 0.1);
    cl[2] = hierarchicalClustering(dtTrain, 0.25);
    cl[3] = hierarchicalClustering(dtTrain, 0.50);
    cl[4] = homogeneousClustering(dtTrain, false);
    cl[5] = homogeneousClustering(dtTrain, true);
    cl[6] = weightedGreedyClustering(dtTrain);
    cl[7] = hierarchicalClustering(dtTrain, (double)cl[6].clusters.size()/(double)dtTrain.I);

    for (int c=0; c<nbClusterings; c++){
      for (int s=0; s<2; s++){
	string modelt = "QOCT";
	if (s==1){
	  modelt += "H";
	}
	for (int selStrat=0; selStrat<4; selStrat++){
	  now = time(0);
	  date = ctime(& now);
	  cout << "Partition " << p << " - " <<typeOfClustering[c]<< " - " << modelt << " - selStrat " << selStrat << " - start : "<< date<<endl;
	
	  try {
	    training_results_cl tr = learningWithClustering(dtTrain, dtValidation, dtTest, cl[c], baseModel::QOCT, s==0, 4,timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1),Nmin,selStrat);
	    file_out.open(filename, ios_base::app);
	    for (int d=0;d<DMAX-1;d++){
	      file_out << datasetName << "," << p << "," << modelt << "," << timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1)<< ",";
	      file_out << typeOfClustering[c] << "," << selStrat<< "," << (double)cl[c].clusters.size()/(double)dtTrain.I << ",";
	      file_out << d+2 << "," << tr.time[d] << "," << tr.nbIter[d] << "," << tr.nbOpti[d] << "," << tr.nbIterCl[d] << ",";
	      file_out << tr.pseudoGap[d] << "," << tr.finalReduction[d] << "," << tr.errTrain[d] << "," << tr.errTest[d] << ",";
	      file_out << tr.bestTree[d] << endl;
	    }
	    file_out.close();
	  } catch(...){
	    file_out.open(filename, ios_base::app);
	    for (int d=0;d<DMAX-1;d++){
	      file_out << datasetName << "," << p << "," << modelt << "," << timeLimitUniv*(int)(s==0)+timeLimitMultiv*(int)(s==1)<< ",";
	      file_out << typeOfClustering[c] << "," << selStrat<< "," << (double)cl[c].clusters.size()/(double)dtTrain.I << ",";
	      file_out << d+2 << ", ERROR" << endl;
	    }
	    file_out.close();
	  }
	}	
      }
    }
  }
}
