// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Tools for calculating Wigner functions from density operators}
#ifndef CPPQEDCORE_QUANTUMDATA_DISTRIBUTIONFUNCTIONS_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_DISTRIBUTIONFUNCTIONS_H_INCLUDED

#include "BlitzArray.h"
#include "ComplexExtensions.h"
#include "MathExtensions.h"
#include "ParsFwd.h"

#include <boost/math/special_functions/factorials.hpp>


namespace quantumdata {


/// Parameter set for scanFunction
struct ParsFunctionScan {
  
  double
    &fLimitXUL, &fLimitYUL, &fLimitXL, &fLimitXU, &fLimitYL, &fLimitYU, &fStep;
    
  int &fCutoff;
  
  ParsFunctionScan(parameters::Table& p, const std::string& mod="");

};


namespace details {

double w(size_t n, double r, size_t k);

}

/// Calculates the Wigner function corresponding to a harmonic-oscillator mode density matrix expressed in Fock basis
/**
 * We use the formula given by Leonhardt \cite leonhardt
 * \f[W[\rho](x,y)=w(r,0)+2\sum_{k=1}^M \real{w(r,k)\,e^{-i\,k\varphi}},\f]
 * with
 * \f[w(r,k)\equiv\sum_{n=0}^{M-k}\left[\frac1\pi(-1)^n\sqrt{\frac{n!}{(n-k)!}}\,e^{-2\,r^2}\lp2\,r\rp^k L_n^k\lp4\,r^2\rp\right]\rho_{n+k,n}.\f]
 * The difference of \f$\sqrt{2}\f$ in front of \f$r\f$ comes from our different definition of the quadratures: \f$x=\frac{a+a^\dagger}2,\;y=\frac{a-a^\dagger}{2i}\f$.
 * The quadratures are recalculated into polar coordinates: \f$x=r\,\cos\lp\varphi\rp\;y=r\,\sin\lp\varphi\rp\f$.
 * 
 * \note The function is wasteful at the moment from several points of view:
 *   - we do not use the recurrence formula (5.114) of the Leonhardt book
 *   - \f$w(r,k)\f$ gets recalculated for each pair of \f$x\;y\f$ leaving unexploited that it depends only on \f$r\f$
 * 
 * \tparam DensityOperatorFunctor a type modelling a unary density operator
 * 
 * \param truncatedDimension if nonzero, this is used as the dimension instead of the actual dimension of `rho`
 * 
 * \todo Should refer `rho` only via some traits class to make the code really generic.
 * 
 */
template<typename DensityOperatorFunctor>
double wignerFunction(const DensityOperatorFunctor& rho, double x, double y, size_t truncatedDimension=0)
{
  using cppqedutils::sqr; using boost::math::factorial;

  struct Helper
  {
    static dcomp _(const DensityOperatorFunctor& rho, size_t dim, double r, size_t k)
    {
      dcomp res(0.);
      for (size_t n=0; n<dim-k; ++n) res+=details::w(n,r,k)*rho(n+k)(n);
      return res;
    }
  };
  
  const double r=sqrt(sqr(x)+sqr(y)), phi=(x || y) ? atan2(y,x) : 0.;

  const size_t dim=truncatedDimension && truncatedDimension<rho.getDimension(0) ? truncatedDimension : rho.getDimension(0);

  double res=real(Helper::_(rho,dim,r,0));
  
  for (size_t k=1; k<dim; ++k) res+=2*real(Helper::_(rho,dim,r,k)*exp(-1i*(k*phi)));

  return res;
}

template<typename DensityOperator>
double qFunction(const DensityOperator& rho, double x, double y, size_t)
{
  using namespace cppqedutils;

  const dcomp alpha(x,y);

  dcomp qComplex;
  
  for (size_t m=0; m<rho.getDimension(); ++m) for (size_t n=0; n<rho.getDimension(); ++n)
    qComplex+=coherentElement(n,alpha)*coherentElement(m,conj(alpha))*rho(m)(n);
  
  return exp(-sqrAbs(alpha))*real(qComplex);
}
  

class WignerFunctionKernelOld
{
public:
  typedef DArray<1> Hermites;
  
  WignerFunctionKernelOld(double x, double y, size_t dim);

  dcomp operator()(size_t m, size_t n) const;
  
private:
  const Hermites hermite_m2x_, hermite_2y_;
  
};



template<typename DensityOperatorFunctor>
double wignerFunctionOld(const DensityOperatorFunctor& rho, double x, double y, size_t truncatedDimension=0)
// NEEDS_WORK should refer rho only via some traits class to make the code really generic.
{
  using namespace cppqedutils; using boost::math::factorial;
  
  const size_t dim=truncatedDimension ? truncatedDimension : rho.getDimension(0);
  
  const WignerFunctionKernelOld wignerFunctionKernel(x,y,dim);
  
  double res=0;

  for (size_t m=0; m<dim; ++m) {
    res+=real(minusOneToThePowerOf(m)/factorial<double>(m)*rho(m,m)*pow(2.*1i,-2*m)*wignerFunctionKernel(m,m));
    
    for (size_t n=m+1; n<dim; ++n)
      res+=2.*real(minusOneToThePowerOf(m)/sqrt(factorial<double>(n)*factorial<double>(m))*rho(n,m)*pow(2.*1i,-m-n)*wignerFunctionKernel(n,m));
    
  }
  
  return 2./PI*exp(-2*(sqr(x)+sqr(y)))*res;
  
}


/// Creates a map of the Wigner function over the region specified by `pfs` and streams it in a format suitable for gnuplot pm3d maps
template<typename DistributionFunctor,typename DensityOperator>
std::ostream& scanFunction(DistributionFunctor distributionFunctor, const DensityOperator& rho, std::ostream& os, const ParsFunctionScan& pfs)
{
  if (pfs.fLimitXUL) pfs.fLimitXL=-(pfs.fLimitXU=pfs.fLimitXUL);
  if (pfs.fLimitYUL) pfs.fLimitYL=-(pfs.fLimitYU=pfs.fLimitYUL);
  for (double x=pfs.fLimitXL; x<pfs.fLimitXU; x+=pfs.fStep) {
    for (double y=pfs.fLimitYL; y<pfs.fLimitYU; y+=pfs.fStep) os<<x<<"\t"<<y<<"\t"<<distributionFunctor(rho,x,y,pfs.fCutoff)<<std::endl;
    os<<std::endl;
  }
  return os;
}



}

#endif // CPPQEDCORE_QUANTUMDATA_DISTRIBUTIONFUNCTIONS_H_INCLUDED
