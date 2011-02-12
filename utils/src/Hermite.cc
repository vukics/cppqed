#include<gsl/gsl_poly.h>

#include "Conversions.h"

#include "Hermite.h"
#include "HermiteCoefficients.h"

namespace cpputils {

using HermiteCoefficients::lim;
using HermiteCoefficients::c  ;

const double Hermite(size_t n, double x)
{

  if (n>=lim) throw HermiteOverflow();

  return gsl_poly_eval(c[n],size2Int(n+1),x);

}

const dcomp C_poly_eval(const double* coeff, const int len, const dcomp& x)
{
  dcomp xtothei=1., result=0.;
  for (int i=0; i<len; (i++, xtothei*=x)) result+=coeff[i]*xtothei;
  return result;
}

const dcomp HermiteC(size_t n, const dcomp& x)
{
  using namespace HermiteCoefficients;

  if (n>=lim) throw HermiteOverflow();

  return C_poly_eval(c[n],size2Int(n+1),x);

}

}
