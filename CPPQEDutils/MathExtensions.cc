// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MathExtensions.h"

#include <boost/math/special_functions/factorials.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include <stdexcept>

  
const double cppqedutils::PI(M_PI);
const double cppqedutils::SQRTPI(M_SQRTPI);
const double cppqedutils::EULER(M_E);

int cppqedutils::sign(double x) {return GSL_SIGN(x);}
int cppqedutils::fcmp(double x, double y, double eps) {return gsl_fcmp(x,y,eps);}

double cppqedutils::sqr(double x) {return gsl_pow_2(x);}

double cppqedutils::sqrAbs(dcomp x) {return sqr(real(x))+sqr(imag(x));} // saves the sqrt

double cppqedutils::fact(unsigned n)
{
  if (n>GSL_SF_FACT_NMAX) throw std::out_of_range("Factorial of"+std::to_string(n));
  return gsl_sf_fact(n);
}

double cppqedutils::choose(unsigned n, unsigned m)
{
  return gsl_sf_choose(n,m);
}

dcomp cppqedutils::coherentElement(unsigned long n, dcomp alpha)
{
  using namespace boost::math;
  return n ? n<max_factorial<double>::value ? pow(alpha,n)/sqrt(factorial<double>(n)) 
                                            : pow(2*n*PI,-.25)*pow(alpha/sqrt(n/EULER),n)
           : 1.;
}
