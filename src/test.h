#pragma once
#include "train_tree.h"
#include "gurobi_c++.h"
#include<iostream>
#include<fstream>
#include<string>
#include <ctime>
using namespace std;

class info_matrix{
 public :
  int nb_columns;
  string* columns_name;
  int nb_lines;
  string* content;

  info_matrix(){}
  info_matrix(int n, int l, string col_names[]);
  void write_line(int line_num, string l[]);
  void write_csv(string namefile, int beg_line=0, int end_l=-1);
};

class statistical_test_matrix : public info_matrix{
 public :
  int nb_models, nb_part, nb_depths;
  string dataset;
  int timeL, Nmin;
  float alpha;

  statistical_test_matrix(){}
  statistical_test_matrix(int nb_mod, int nb_p, int nb_dep, string dts, int Nm, float alph, vector<double> tl);
  void write_line(int p, int d, int model_cpt, string model, solution sol, string error_test, string error_test_pp);
  void write_result_partition(int p);
};

class learning_test_matrix : public info_matrix{
 public :
  int nb_models, nb_part, DMAX;
  string dataset;

  learning_test_matrix(){}
  learning_test_matrix(int nb_mod, int nb_p, string dts, int D);
  void write_line(int p, int model_cpt, string model, double time_l, training_results tr);
  void write_result_partition(int p);
};

void test(string dataset);
void learning_test(string dataset, int depth_max);
void testClust(string datasetName, float time_l);
void clustStats(string datasetName);
