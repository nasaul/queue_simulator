// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_randist.h>
#include "FunctionsRandom.h"
#include "RandomNumbers.h"

//--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++
double RanMT(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
}
//--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++
double ExpoMT(const double &beta)
{
    double uni;
    uni = RanMT();
    return (-beta * log(uni));

}
//--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++
int multinom(arma::mat trans_matrix, int state, const gsl_rng * r, int K){
  double prob[K];
  for(int i = 0; i < K; i++){
    prob[i] = trans_matrix(state - 1, i);
  }
  unsigned int x[K];
  gsl_ran_multinomial(r, K, 1, prob, x);
  int val = 0;
  for(int i = 0; val == 0; i++){
    if(x[i] == 1){
      val = i + 1;
    }
  }
  return(val);
}
