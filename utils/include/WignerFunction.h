// -*- C++ -*-
#ifndef UTILS_INCLUDE_WIGNERFUNCTION_H_INCLUDED
#define UTILS_INCLUDE_WIGNERFUNCTION_H_INCLUDED

#include "BlitzArray.h"
#include "ComplexExtensions.h"
#include "MathExtensions.h"


namespace cpputils {


/*
class BufferedHermitePolynomial
{
public:
  typedef std::vector<double> Values;
  
  typedef std::map<double,Values,boost::function<bool(double,double)> > Impl;
  // use approximate comparison gsl_fcmp
  
  
  // @ construction, read Impl from file => query and append it as needed => @ destruction, replace IF changed
  // use serialization

  long double operator()(unsigned n, double x) const;
  
private:
  mutable Impl;
};
*/
  


class WignerFunctionKernel
{
public:
  typedef TTD_DARRAY(1) Hermites;
  
  WignerFunctionKernel(double x, double y, size_t dim);

  dcomp operator()(size_t m, size_t n) const;
  
private:
  const Hermites hermite_2x_, hermite_m2y_;
  
};



template<typename DensityOperator>
double wignerFunction(const DensityOperator& rho, double x, double y, size_t truncatedDimension=0)
// NEEDS_WORK should refer rho only via some traits class to make the code really generic.
{
  using namespace std; using namespace mathutils;
  
  const size_t dim=truncatedDimension ? truncatedDimension : rho.getDimensions(0);
  
  const WignerFunctionKernel wignerFunctionKernel(x,y,dim);
  
  double res=0;

  for (int m=0; m<dim; ++m) {
    res+=real(rho(m,m)/fact(m)*pow(-1.,m)*pow(2.*DCOMP_I,-2*m)*wignerFunctionKernel(m,m));
    
    for (int n=m+1; n<dim; ++n)
      res+=2.*real(rho(n,m)/sqrt(fact(n)*fact(m))*pow(-1.,m)*pow(2.*DCOMP_I,-m-n)*wignerFunctionKernel(n,m));
    
  }
  
  return 2./PI*exp(-2*(sqr(x)+sqr(y)))*res;
  
}


/*
template<typename DensityOperator, typename Kernel>
double wignerFunction(const DensityOperator& rho, const dcomp& alpha, const Kernel&);

+ a bunch of further overloads

*/

}

#endif // UTILS_INCLUDE_WIGNERFUNCTION_H_INCLUDED
