// -*- C++ -*-
#ifndef UTILS_INCLUDE_INTEGRATION_H_INCLUDED
#define UTILS_INCLUDE_INTEGRATION_H_INCLUDED

#include "ComplexExtensions.h"

namespace cpputils {

double RealIntegration(double (*)(double, void*), double, double, void*, double epsabs=1e-6, double epsrel=1e-3);

}

#endif // UTILS_INCLUDE_INTEGRATION_H_INCLUDED
