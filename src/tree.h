#pragma once
#include "data_formatting.h"
#include "gurobi_c++.h"

class Tree
{    
public:
    // paramètres de l'abre
    int D, J, K;
    //vector<int>* A_L;
    //vector<int>* A_R;
    vector<int>* leaves;
    int* right_most_leaf; // RML[t] forall t \in NUL, RML[t]=t forall t \in L et RML[t] = right most leaf forall t \in N 
    int L,N,E;
    // définition des règles de branchements
    double* a;
    double* b;
    double* eps;
    // attribution des classes aux feuilles
    int* c; // les classes seront entre 0 et K-1

    // fonctions de création d'un arbre
    Tree();
    Tree(int D_, int J_, int K_); // initialisation sans écriture des valeurs de a, b etc
    
    Tree copy();
    // fonctions utilitaire pour la construction d'un arbre
    void get_tree_structure();

    void correctSplitFunctions();  // correcting split functions when coefficent are too small
    // fonction permettant de "remplir" les coefficients d'un arbre défini par Tree(int D_, int J_, int K_)
    void post_processing_b(dataset& dt, bool missClassifCounts=false, float rho=0.5); // remplir les valeurs de a, b etc avec le post-traitement sur b
    void post_processing_a_b(dataset& dt, bool missClassifCounts=false); // remplir les valeurs de a, b etc avec le post-traitement sur a et b
    void compute_ILIR(vector<int>* I_L, vector<int>* I_R, dataset& dt, bool missClassifCounts=false); // IL and IR must be empty

    // fonction permettant d'agrandir un arbre
    Tree bigger_tree(int new_D);
    
    // fonction qui pousse b à son max
    Tree adjustTree(dataset& dt);

    // fonctions qui permettent d'enlever un coeff à l'arbre
    Tree reduceComplexity(dataset& dt, int Nmin=0, double mu = 0.0001); // en resolvant un mini-programme
    bool changeSplit(int t, vector<int> points, int& misclassif, dataset& dt, int Nmin=0, double mu=0.0001);
    Tree removeSplit(dataset& dt); // en enlevant un split

    int countComplexity();
    
    // fonction permettant de prédire la classe de données
    int predict_class(double x[]);
    int predict_leaf(double x[]);
    void predict_classes(dataset& dt,int predictions[]);
    int prediction_errors(dataset& dt);
    void data_points_per_leaves(dataset& dt, int repartition[]); // repartition[t*K+k] is the number of data point of class k in node t
    void data_points_in_last_split(dataset& dt, vector<int> points_in_node[]);
    void predict_leaves(dataset& dt, int leaves[]);
    void predict_paths(dataset& dt, vector<int> paths[]); // for F warm-start

    // fonction permettant d'écrire ou lire un arbre avec des fichiers .txt
    void write_tree(string namefile);
    Tree(string namefile); // équivalent de la fonction read_tree
};
