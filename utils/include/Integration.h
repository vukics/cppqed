// -*- C++ -*-
#ifndef _INTEGRATION_H
#define _INTEGRATION_H

#include "ComplexExtensions.h"

namespace cpputils {

double RealIntegration(double (*)(double, void*), double, double, void*, double epsabs=1e-6, double epsrel=1e-3);

}

#endif // _INTEGRATION_H
