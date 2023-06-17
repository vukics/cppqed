// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Mode.h"


mode::MultiDiagonal mode::aOp(size_t cutoff)
{
  MultiDiagonal res{};
  res.diagonals.insert({MultiDiagonal::Offsets{1},MultiDiagonal::Diagonal{{cutoff-1}, [c=cutoff] (size_t) {
    // TODO: pending support for std::ranges::to, the following solution will be possible:
    // return std::ranges::to<MultiDiagonal::Diagonal::StorageType>(std::views::iota(0uz, c-1) | std::views::transform([] (size_t num) -> dcomp { return std::sqrt(num+1); }));
    auto res{::quantumdata::noInit<1>(c-1)};
    for (size_t i=0; i<res.size(); ++i) res[i]=sqrt(double(i+1));
    return res;
  }}});
  return res;
}


mode::MultiDiagonal mode::aDagOp(size_t cutoff)
{
  auto res(aOp(cutoff)); res.hermitianConjugateSelf(); return res;
}


mode::MultiDiagonal mode::nOp(size_t cutoff)
{
  return compose(aDagOp(cutoff),aOp(cutoff));
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
  for (size_t n=0; n<psi.ubound(0); ++n) psi(n)=fact*sqrt(n+1)*psi(n+1);
  psi(psi.ubound(0))=0;
}
    

void mode::aDagJump(StateVectorLow& psi, double fact)
{
  for (size_t n=psi.ubound(0); n>0; --n) psi(n)=fact*sqrt(n)*psi(n-1);
  psi(0)=0;
}


void mode::aSuperoperator(const DensityOperatorLow& rho, DensityOperatorLow& drhodt, double fact)
{
  for (int m=0; m<rho.ubound(0); ++m)
    for (int n=0; n<rho.ubound(1); ++n)
      drhodt(m,n)+=fact*sqrt((m+1)*(n+1))*rho(m+1,n+1);
}


void mode::aDagSuperoperator(const DensityOperatorLow& rho, DensityOperatorLow& drhodt, double fact)
{
  for (int m=1; m<=rho.ubound(0); ++m)
    for (int n=1; n<=rho.ubound(1); ++n)
      drhodt(m,n)+=fact*sqrt(m*n)*rho(m-1,n-1);
}


Lindblad<1> mode::photonLoss(double kappaTimes_nThPlus1)
{
  return {
    .label{"excitation loss"},
    .jump{[=](StateVectorLow& psi) {aJump(psi,sqrt(2.*kappaTimes_nThPlus1));}},
    .rate{[=](const LazyDensityOperator& m) {return 2.*kappaTimes_nThPlus1*photonnumber(m);}},
    .superoperator{[=](const DensityOperatorLow& rho, DensityOperatorLow& drhodt) {aSuperoperator(rho,drhodt,2.*kappaTimes_nThPlus1);}}
  };
}

Lindblad<1> mode::photonGain(double kappaTimes_nTh)
{
  return {
    .label{"excitation absorption"},
    .jump{[=](StateVectorLow& psi) {aDagJump(psi,sqrt(2.*kappaTimes_nTh));}},
    .rate{[=](const LazyDensityOperator& m) {return 2.*kappaTimes_nTh*(photonnumber(m)+m.trace());}},
    .superoperator{[=](const DensityOperatorLow& rho, DensityOperatorLow& drhodt) {aDagSuperoperator(rho,drhodt,2.*kappaTimes_nTh);}}
  };
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
