#include "utilities.h"

bool lessThan(double a, double b){
  return a+GLOBAL_PRECISION<b;
}

bool lessThanOrEqualTo(double a, double b){
  return a-GLOBAL_PRECISION<b;
}

bool equalTo(double a, double b){
  return (b-a < GLOBAL_PRECISION) and (a-b < GLOBAL_PRECISION);
}

bool greaterThanOrEqualTo(double a, double b){
  return a + GLOBAL_PRECISION > b;
}

bool greaterThan(double a, double b){
  return a - GLOBAL_PRECISION > b;
}
