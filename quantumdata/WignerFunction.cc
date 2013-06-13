#include "WignerFunction.h"

#include "Pars.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/laguerre.hpp>


using namespace mathutils;
using namespace boost::math;


namespace quantumdata {


ParsWignerFunctionScan::ParsWignerFunctionScan(parameters::ParameterTable& p, const std::string& mod)
  : wfLimitXUL(p.addTitle("Wigner function scan",mod).add("wfLimitXUL","",2.)),
    wfLimitYUL(p.add("wfLimitYUL","",2.)),
    wfLimitXL(p.add("wfLimitXL","",-2.)),
    wfLimitXU(p.add("wfLimitXU","",2.)),
    wfLimitYL(p.add("wfLimitYL","",-2.)),
    wfLimitYU(p.add("wfLimitYU","",2.)),
    wfStep(p.add("wfStep","",.1)),
    wfCutoff(p.add("wfCutoff","",100))
{}


  
namespace {


const WignerFunctionKernelOld::Hermites fillWithHermite(size_t dim, double x)
{
  WignerFunctionKernelOld::Hermites res(2*dim-1);
  res(0)=hermite(0,x); res(1)=hermite(1,x);
  for (unsigned l=1; l<res.size()-1; ++l)
    res(l+1)=hermite_next(l,x,res(l),res(l-1));
  return res;
}


}

double details::w (size_t n, double r, size_t k)
{
  const double sqrR=sqr(r);
  return minusOneToThePowerOf(n)/PI*sqrt(factorial<double>(n)/factorial<double>(n+k))*exp(-2*sqrR)*pow(2*r,k)*laguerre(n,k,4*sqrR);
}


WignerFunctionKernelOld::WignerFunctionKernelOld(double x, double y, size_t dim)
  : hermite_m2x_(fillWithHermite(dim,-2*x)), hermite_2y_(fillWithHermite(dim,2*y))
{}


dcomp WignerFunctionKernelOld::operator()(size_t m, size_t n) const
{
  dcomp res(0);

  for (size_t u=0; u<=m; ++u) for (size_t v=0; v<=n; ++v)
    res+=
      binomial_coefficient<double>(m,u)*
      binomial_coefficient<double>(n,v)*
      minusOneToThePowerOf(v)*
      pow(DCOMP_I,u+v)*
      hermite_m2x_(    u+v)*
      hermite_2y_ (n+m-u-v);

  return res;

}


} // quantumdata

