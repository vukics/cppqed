#include "FuzzyDouble.h"

#include<gsl/gsl_math.h>

#include<limits>


namespace cpputils {


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


} // cpputils
