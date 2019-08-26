#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "FunctionsRandom.h"
#include "SupportFunctions.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
/* ----------------------------------------------------------------
 Crea matriz de transición para cadena de markov
 ----------------------------------------------------------------*/
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
/* ----------------------------------------------------------------
 Genera cadena de Markov
 ----------------------------------------------------------------*/
// [[Rcpp::export]]
Rcpp::IntegerVector sim_mc(int P, arma::mat transition_matrix, int seed, int init){
  int states         = transition_matrix.n_cols;
  gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, seed);
  Rcpp::IntegerVector sim(P);
  sim(0) = init;
  for(int i = 1; i < P; i++){
    sim(i) = multinom(transition_matrix, sim(i - 1), rand_gen, states);
  }
  gsl_rng_free(rand_gen);
  return(sim);
}
