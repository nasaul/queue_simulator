#ifndef FunctionsFunctionsRandom_h
#define FunctionsFunctionsRandom_h

#include <stdio.h>
#include <Rcpp.h>

double RanMT(void);
double ExpoMT(Rcpp::NumericVector params);
// double UnifMT(double &a, double &b);
// double NormMT(double &media,double &des);
#endif /* Functions_h */
