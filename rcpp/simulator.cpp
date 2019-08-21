// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ostream>
#include <iostream>
using namespace std;

// [[Rcpp::export]]
arma::mat create_matrix (int n){
  arma::Mat<double> trans_matrix = arma::mat(n, n, arma::fill::randu);
  arma::colvec      row_sum = arma::sum(trans_matrix, 1);
  arma::Mat<double> result = arma::mat(n, n, arma::fill::zeros);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      result(i, j) =  trans_matrix(i, j) / row_sum(i);
    }
  }
  return(result);
}

int multinom(arma::mat trans_matrix, int state, const gsl_rng * r, int K){
  double prob[K];
  for(int i = 0; i < K; i++){
    prob[i] = trans_matrix(state - 1, i);
  }
  Rcpp::IntegerVector x(K);
  gsl_ran_multinomial(r, K, 1, prob, (unsigned int *) x.begin());
  int val;
  for(int i = 0; i < K; i++){
    if(x(i) == 1){
      val = i + 1;
    }
  }
  return(val);
}

// [[Rcpp::export]]
Rcpp::IntegerVector sim_mc(int P, arma::mat transition_matrix, int seed){
  int states = transition_matrix.n_cols;
  Rcpp::IntegerVector sim(P);
  gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, seed);
  sim(0) = 1;
  for(int i = 1; i < P; i++){
    sim(i) = multinom(transition_matrix, sim(i - 1), rand_gen, states);
  }
  gsl_rng_free(rand_gen);
  return(sim);
}
