// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;
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