//
//  Functions.hpp
//  TestFunctionsRcpp
//
//  Created by Luis Moncayo on 8/23/19.
//  Copyright © 2019 Luis Moncayo. All rights reserved.
//

#ifndef FunctionsFunctionsRandom_h
#define FunctionsFunctionsRandom_h

#include <stdio.h>

double RanMT(void);
double ExpoMT(const double &beta);
bool ExpoVec(long &NExpo, double Tiempos[], double &Media);
int multinom(arma::mat trans_matrix, int state, const gsl_rng * r, int K);

#endif /* Functions_h */
