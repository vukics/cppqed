// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MathExtensions.h"

#include <boost/math/special_functions/factorials.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

namespace mathutils {
  
const double PI(M_PI);
const double SQRTPI(M_SQRTPI);
const double EULER(M_E);

int sign(double x) {return GSL_SIGN(x);}
int fcmp(double x, double y, double eps) {return gsl_fcmp(x,y,eps);}

double sqr(double x) {return gsl_pow_2(x);}

double sqrAbs(const dcomp& x) {return sqr(real(x))+sqr(imag(x));} // saves the sqrt

double fact(unsigned n) throw(FactOverflow)
{
  if (n>GSL_SF_FACT_NMAX) throw FactOverflow();
  return gsl_sf_fact(n);
}

double choose(unsigned n, unsigned m)
{
  return gsl_sf_choose(n,m);
}

dcomp coherentElement(unsigned long n, const dcomp& alpha)
{
  using namespace boost::math;
  return n ? n<max_factorial<double>::value ? pow(alpha,n)/sqrt(factorial<double>(n)) 
                                            : pow(2*n*PI,-.25)*pow(alpha/sqrt(n/EULER),n)
           : 1.;
}


} // mathutils
