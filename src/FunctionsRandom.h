#ifndef FunctionsFunctionsRandom_h
#define FunctionsFunctionsRandom_h

#include <stdio.h>
#include <Rcpp.h>

double RanMT(void);
double ExpoMT(Rcpp::NumericVector params);
double UnifMT(Rcpp::NumericVector params);
double NormMT(Rcpp::NumericVector params);
#endif /* Functions_h */
