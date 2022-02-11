// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "FuzzyDouble.h"

#include<gsl/gsl_math.h>

#include<limits>


namespace cppqedutils {


const double FuzzyDouble::eps=10.*std::numeric_limits<double>::epsilon();


bool operator==(FuzzyDouble a, FuzzyDouble b)
{
  return !gsl_fcmp(a,b,FuzzyDouble::eps);
}


bool operator< (FuzzyDouble a, FuzzyDouble b)
{
  return gsl_fcmp(a,b,FuzzyDouble::eps)<0;
}


bool operator==(FuzzyDouble a, double b)
{
  return !gsl_fcmp(a,b,FuzzyDouble::eps);
}


bool operator< (FuzzyDouble a, double b)
{
  return gsl_fcmp(a,b,FuzzyDouble::eps)<0;
}


bool operator> (FuzzyDouble a, double b)
{
  return gsl_fcmp(a,b,FuzzyDouble::eps)>0;
}


} // cppqedutils
