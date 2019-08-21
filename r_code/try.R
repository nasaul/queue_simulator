library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(purrr)

Rcpp::sourceCpp("rcpp/simulator.cpp")
a = create_matrix(100)
sim_mc(1000, a, 7, 8)

create_matrix_r <- function(n){
  mat  <-  matrix(ncol = n, nrow = n)
  mat_ind <- map(
    1:n,
    function(x){
      x <- runif(n)
      y <- x / sum(x)
    }
  )
  for(i in 1:n){
    mat[i, ] = mat_ind[[i]]
  }
  gc()
  return(mat)
}

bench::mark(
  create_matrix(1000),
  create_matrix_r(1000),
  check = FALSE,
  iterations = 50
  # filter_gc = F
) %>% View
