// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Mode.h"

#include <boost/math/special_functions/factorials.hpp>


using namespace structure;


mode::MultiDiagonal mode::aOp(size_t cutoff)
{
  MultiDiagonal res;
  res.diagonals[1].emplace(MultiDiagonal::Offsets{1},MultiDiagonal::Diagonal{ [c=cutoff] {
    // TODO: pending support of std::ranges::to, the following solution will be possible:
    // return std::ranges::to<MultiDiagonal::Diagonal::StorageType>(std::views::iota(0uz, c-1) | std::views::transform([] (size_t num) -> dcomp { return std::sqrt(num+1); }));
    auto res{noInit(c-1)};
    for (size_t i=0; i<res.size(); ++i) res[i]=sqrt(double(i+1));
    return res;
  } () } );
  return res;
}


mode::MultiDiagonal mode::aDagOp(size_t cutoff)
{
  return hermitianConjugateOf(aOp(cutoff));
}


mode::MultiDiagonal mode::nOp(size_t cutoff)
{
  return aDagOp(cutoff) | aOp(cutoff);
}


TimeIndependentJump<1> mode::aJump(double fact)
{
  return [=] (StateVectorView<1> psi) {
    for (size_t n=0; n<psi.extents[0]-1; ++n) psi(n)=fact*sqrt(n+1)*psi(n+1);
    psi(psi.extents[0]-1)=0;
  };
}


TimeIndependentJump<1> mode::aDagJump(double fact)
{
  return [=] (StateVectorView<1> psi) {
    for (size_t n=psi.extents[0]-1; n>0; --n) psi(n)=fact*sqrt(n)*psi(n-1);
    psi(0)=0;
  };
}


TimeIndependentSuperoperator<1> mode::aSuperoperator(double fact)
{
  return [=] (DensityOperatorConstView<1> rho, DensityOperatorView<1> drhodt) {
    for (int m=0; m<rho.extents[0]-1; ++m)
      for (int n=0; n<rho.extents[1]-1; ++n)
        drhodt(m,n)+=fact*sqrt((m+1)*(n+1))*rho(m+1,n+1);
  };
}


TimeIndependentSuperoperator<1> mode::aDagSuperoperator(double fact)
{
  return [=] (DensityOperatorConstView<1> rho, DensityOperatorView<1> drhodt) {
    for (int m=1; m<=rho.extents[0]-1; ++m)
      for (int n=1; n<=rho.extents[1]-1; ++n)
        drhodt(m,n)+=fact*sqrt(m*n)*rho(m-1,n-1);
  };
}


Lindblad<1> mode::photonLoss(double kappa, double nTh)
{
  double fact=2.*kappa*(nTh+1);
  return {
    .label{"photon loss"},
    .jump{aJump(sqrt(fact))},
    .rate{ [=] (StateVectorConstView<1> psi) {return fact*photonnumber(psi);} },
    .superoperator{aSuperoperator(sqrt(fact))}
  };
}


Lindblad<1> mode::photonGain(double kappa, double nTh)
{
  double fact=2.*kappa*nTh;
  return {
    .label{"photon gain"},
    .jump{aDagJump(sqrt(fact))},
    .rate{ [=] (StateVectorConstView<1> psi) {
      double res=0;
      for (size_t n=0; n<psi.extents[0]; ++n) res+=(n+1)*sqrAbs(psi(n));
      return fact*res;
    }},
    .superoperator{aDagSuperoperator(sqrt(fact))}
  };
}


dcomp coherentElement(unsigned long n, dcomp alpha)
{
  using namespace boost::math;
  return n ? n<max_factorial<double>::value ? pow(alpha,n)/sqrt(factorial<double>(n))
                                            : pow(2*n*std::numbers::pi,-.25)*pow(alpha/sqrt(n/std::numbers::e),n)
           : 1.;
}

StateVector<1> mode::coherent(dcomp alpha, size_t dim)
{
  return { {dim} , [=] (size_t e) {
    auto r{noInit(e)};
    double norm=exp(-sqrAbs(alpha)/2.);
    for (size_t n=0; n<dim; ++n) r[n]=norm*coherentElement(n,alpha);
    return r;
  } };
}


StateVector<1> mode::fock(size_t n, size_t dim, double phase)
{
  if (n>=dim) throw std::overflow_error("Fock state "+std::to_string(n)+" higher than dim "+std::to_string(dim));
  StateVector<1> res({dim},zeroInit);
  res(n)=exp(1i*phase);
  return res;
}


StateVector<1> mode::init(const Pars& p)
{
  return p.initFock ? fock(p.initFock,p.cutoff) : coherent(p.init,p.cutoff);
}
