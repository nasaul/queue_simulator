library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)

Rcpp::sourceCpp("rcpp/simulator.cpp")
a = create_matrix(100)
sim_mc(10000, a, 1)
