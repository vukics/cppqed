// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Mode.h"


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


/*
double mode::photonnumber(const LazyDensityOperator& matrix)
{
  double res{0};
  for (size_t n=0; n<matrix.getDimension(); ++n) res+=n*matrix(n);
  return res;
}


void mode::aJump(StateVectorLow& psi, double fact)
{
}
    

void mode::aDagJump(StateVectorLow& psi, double fact)
{
}




ExpectationValue<1> mode::photonnumberEV_Variance{
  .label{"<number operator>","VAR(number operator)"},
  .process{calculateVariance},
  .eva{[] (const LazyDensityOperator& m) {
    EV_Array res{0.,2};
    for (int n=1; n<int(m.getDimension()); n++) {
      res[0]+=  n*m(n);
      res[1]+=n*n*m(n);
    }    
    return res;
}}};


ExpectationValue<1> mode::ladderOperatorEV {
  .label{"real(<ladder operator>)","imag(\")"},
  .eva{[] (const LazyDensityOperator& m) {
    EV_Array res{0.,2};
    for (int n=1; n<int(m.getDimension()); n++) {
      dcomp offdiag(sqrt(n)*m(n)(n-1));
      res[0]+=real(offdiag);
      res[1]+=imag(offdiag);
    }    
    return res;
}}};


ExpectationValue<1> mode::monitorCutoff {
  .label{"|Psi(cutoff-1)|^2"},
  .eva{[] (const LazyDensityOperator& m) {
    EV_Array res{0.,1};
    res[0]=m(m.getDimension()-1);
    return res;
}}};




const Tridiagonal aop(size_t dim)
{
  typedef Tridiagonal::Diagonal Diagonal;
  Diagonal diagonal(dim-1);
  return Tridiagonal(Diagonal(),1,Diagonal(),diagonal=sqrt(blitz::tensor::i+1.));
}

*/
