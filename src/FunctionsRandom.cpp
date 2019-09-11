#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
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
