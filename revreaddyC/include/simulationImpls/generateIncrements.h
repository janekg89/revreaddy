#ifndef INCREMENTS_INCLUDED
#define INCREMENTS_INCLUDED

#include <boost/multi_array.hpp>
#include <complex>
/*#include </home/mi/janekg89/Dokumente/Masterarbeit/test_cython/include/fftw3.h>*/
#include <fftw3.h>
#include<stdio.h>
#include <iostream>
#include<stdlib.h>
#include <algorithm>
#include "Random.h"

namespace janek {
    boost::multi_array<double,2> generateIncrements(unsigned long N, double& D,double& tau, double& alpha, Random * random);
}

#endif // INCREMENTS_INCLUDED

