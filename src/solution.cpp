#include "solution.h"

solution::solution(Tree T_, double o, int et, int nb, float t, double g, int n, double r){
    T =T_;
    obj = o;
    error_train = et;
    nb_br = nb;
    time = t;
    gap = g;
    nodes = n;
    root_rel = r;
}

bool dominating_trees::add_tree(optimal_tree& opt){

  if (nb_sol == 0){
    trees[0] = opt;
    alph[0] = 0.0;
    alph[1] = Lh;
    nb_sol = 1;
    return true;
  }

  int before = -1,
    after = -1,
    dom=0;

  for (int t=0; t<nb_sol; t++){
    double diff_en_alph_min = trees[t].E - opt.E + alph[t]*(trees[t].C - opt.C),
      diff_en_alph_max  = trees[t].E - opt.E + alph[t+1]*(trees[t].C - opt.C);

    if (diff_en_alph_min != diff_en_alph_max){ // pas parallèle
      if (diff_en_alph_min ==0){// elles se coupent en alpha_min
	if (diff_en_alph_max >0){
	  before = t-1;
	}
      }
      else{
	if (diff_en_alph_min >0){
	  if (diff_en_alph_max ==0){
	    after = t+1;
	  }
	  else{
	    if (diff_en_alph_max <0){
	      after = t;
	    }
	    else{ // diff_en_alph_max >0
	      dom += 1;
	    }
	  }
	}
	else{ // diff_en_alph_min <0
	  if (diff_en_alph_max >0){
	    before = t;
	  }
	}
      }
 
    }
    else{// si elles sont parallèles
      if (diff_en_alph_max >0){// elle est en dessous
	dom += 1;
      }
      else{ // elle est au dessus (dominée)
	break;
      }
    }
  }

  if (before == -1){
    if (after == -1){// la solution domine toutes les autres ou aucune
      if (dom==0){ // aucune solution n'est dominée
	return false;
      }
      else{
	trees[0] = opt;
	nb_sol = 1;
	alph[0] = 0.0;
	alph[1] = Lh;
      }
    }
    else if (after ==0){
      double meeting_point = (double)(opt.E - trees[after].E)/(double)(trees[after].C - opt.C);
      for (int t=0; t<nb_sol; t++){
	trees[nb_sol-t] = trees[nb_sol-t-1];
	alph[nb_sol-t+1] = alph[nb_sol-t];
      }
      alph[1] = meeting_point;
      trees[0] = opt;
      nb_sol += 1;
    }
    else{
      double meeting_point = (double)(opt.E - trees[after].E)/(double)(trees[after].C - opt.C);
      alph[1] = meeting_point;
      for (int t=after; t<nb_sol; t++){
	trees[t-after+1] = trees[t];
	alph[t-after+2] = alph[t+1];
      }
      nb_sol = nb_sol - after + 1;
      trees[0] = opt;
    }
  }
  else{// before >= 0
    double meeting_point1 = (double)(opt.E - trees[before].E)/(double)(trees[before].C - opt.C);
    if (after == -1){
      alph[before+1] = meeting_point1;
      alph[before+2] = Lh;
      trees[before+1] = opt;
      nb_sol = before + 2;
    }
    else{
      double meeting_point2 = (double)(opt.E - trees[after].E)/(double)(trees[after].C - opt.C);
      if (dom > 0){ // ie on va supprimer des sol
	trees[before+1] = opt;
	alph[before+1] = meeting_point1; 
	alph[before+2] = meeting_point2;
	
	for (int t=after; t<nb_sol; t++){
	  trees[before + 2 + t-after] = trees[t];
	  alph[before + 2 + t-after+1] = alph[t+1];
	}
	nb_sol = before + 2 + nb_sol - after;
      }
      else{ // on ajoute une solution
	for (int t=0; t< nb_sol-after;t++){
	  trees[nb_sol-t] = trees[nb_sol-t-1];
	  alph[nb_sol-t+1] = alph[nb_sol-t];
	}
	trees[before+1] = opt;
	alph[before+1] = meeting_point1;
	alph[before+2] = meeting_point2;
	nb_sol += 1;
      }
    }
  }
  return true;
}

int dominating_trees::BestTree(bool smallest_alpha){
  int best = 0,
    score = trees[0].Ev;
  for (int i=1; i<nb_sol; i++){
    bool new_best = false;
    if (smallest_alpha){
      new_best = trees[i].Ev<score;
    }
    else{
      new_best = trees[i].Ev<=score;
    }
    if (new_best){
      best=i;
      score = trees[i].Ev;
    }
  }
  return best;
}

int dominating_trees::BestWarmstart(int D, int C){
  int best = -1,
    score = Lh;
  for (int i=0; i<nb_sol; i++){
    if (trees[i].D <= D && trees[i].C <= C){
      if (trees[i].E < score){
	best = i;
	score = trees[i].E;
      }
    }
  }
  return best;
}
