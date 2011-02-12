#include "Integration.h"

#include<gsl/gsl_integration.h>

double cpputils::RealIntegration(double (*f)(double, void*), double xl, double xu, void* params, double epsabs, double epsrel)
{
  const gsl_function ff={f,params};
  double result, abserr;
  size_t neval;

  gsl_integration_qng(&ff,xl,xu,epsabs,epsrel,&result,&abserr,&neval);

  return result;

}


