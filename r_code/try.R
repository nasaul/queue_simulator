run.mc.sim <- function( P, num.iters = 50 ) {

  # number of possible states
  num.states <- nrow(P)

  # stores the states X_t through time
  states     <- numeric(num.iters)
  print(states)
  # initialize variable for first state
  states[1]    <- 1
  # print(states)
  for(t in 2:num.iters) {

    # probability vector to simulate next state X_{t+1}
    p  <- P[states[t-1], ]
    states[t] <-  which(rmultinom(1, 1, p) == 1)
  }
  return(states)
}

Rcpp::sourceCpp("rcpp/simulator.cpp")
a = create_matrix(10000)
run.mc.sim(create_matrix(50))
dim(a)
