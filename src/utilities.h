#include "GurobiEnvironment.h"

const double mu = 0.0001; // 10-4
const double GLOBAL_PRECISION = 0.000001; // 10-6
const bool setPrecision = true;

/*
const double mu = 0.00001;
const double GLOBAL_PRECISION = mu/10;
*/

bool lessThan(double a, double b);
bool lessThanOrEqualTo(double a, double b);
bool equalTo(double a, double b);
bool greaterThanOrEqualTo(double a, double b);
bool greaterThan(double a, double b);
