// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "DistributionFunctions.h"

#include "Pars.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/factorials.hpp>


using namespace mathutils;
using namespace boost::math;


namespace quantumdata {


ParsFunctionScan::ParsFunctionScan(parameters::ParameterTable& p, const std::string& mod)
  : fLimitXUL(p.addTitle("Distribution function scan",mod).add("fLimitXUL","",2.)),
    fLimitYUL(p.add("fLimitYUL","",2.)),
    fLimitXL(p.add("fLimitXL","",-2.)),
    fLimitXU(p.add("fLimitXU","",2.)),
    fLimitYL(p.add("fLimitYL","",-2.)),
    fLimitYU(p.add("fLimitYU","",2.)),
    fStep(p.add("fStep","",.1)),
    fCutoff(p.add("fCutoff","",100))
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

double details::w(size_t n, double r, size_t k)
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

