#include "WignerFunction.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/hermite.hpp>


using namespace boost::math;


namespace {

const cpputils::WignerFunctionKernel::Hermites fillWithHermite(size_t dim, double x)
{
  cpputils::WignerFunctionKernel::Hermites res(2*dim-1);
  res(0)=hermite(0,x); res(1)=hermite(1,x);
  for (unsigned l=1; l<res.size()-1; ++l)
    res(l+1)=hermite_next(l,x,res(l),res(l-1));
  return res;
}

}


cpputils::WignerFunctionKernel::WignerFunctionKernel(double x, double y, size_t dim)
  : hermite_2x_(fillWithHermite(dim,2*x)), hermite_m2y_(fillWithHermite(dim,-2*y))
{}


dcomp cpputils::WignerFunctionKernel::operator()(size_t m, size_t n) const
{
  using namespace mathutils;

  dcomp res(0);

  for (size_t u=0; u<=m; ++u) for (size_t v=0; v<=n; ++v)
    res+=
      binomial_coefficient<double>(m,u)*
      binomial_coefficient<double>(n,v)*
      pow(-1.,int(v))*
      pow(DCOMP_I,u+v)*
      hermite_2x_ (n+m-u-v)*
      hermite_m2y_(    u+v);

  return res;

}