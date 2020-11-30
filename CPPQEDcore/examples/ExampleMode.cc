// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ExampleMode.h"

#include "ElementLiouvillean.h"
#include "Tridiagonal.tcc"
#include "TridiagonalHamiltonian.tcc"

#include <boost/bind.hpp>

using boost::bind;


void aJump   (StateVectorLow&, double kappa_nPlus1);
void aDagJump(StateVectorLow&, double kappa_n     );

double aJumpRate   (const LazyDensityOperator&, double kappa_nPlus1);
double aDagJumpRate(const LazyDensityOperator&, double kappa_n     );


basic::PumpedLossyMode::PumpedLossyMode(double delta, double kappa, dcomp eta, double nTh, size_t cutoff)
  : Free(cutoff,{},
         {
           CF{"(kappa*(2*nTh+1),delta)",dcomp(kappa*(2*nTh+1),delta),     cutoff },
           CF{"eta"                    ,eta                         ,sqrt(cutoff)}
         }),
    TridiagonalHamiltonian<1,false>(dcomp(-kappa*(2*nTh+1),delta)*nop(cutoff)
                                    +
                                    tridiagPlusHC_overI(conj(eta)*aop(cutoff))),
    ElementLiouvilleanStrategies<1,2>(JumpStrategies(bind(aJump   ,_1,kappa*(nTh+1)),
                                                     bind(aDagJump,_1,kappa* nTh   )),
                                      JumpRateStrategies(bind(aJumpRate   ,_1,kappa*(nTh+1)),
                                                         bind(aDagJumpRate,_1,kappa* nTh   )),
                                      "Mode",{"photon loss","photon absorption"}),
    ElementAveraged<1>("PumpedLossyMode",{"<number operator>","real(<ladder operator>)","imag(\")"})
{
  getParsStream()<<"Pumped lossy mode";
}


auto basic::PumpedLossyMode::average_v(NoTime, const LazyDensityOperator& matrix) const -> const Averages
{
  auto averages(initializedAverages());

  averages(0)=aJumpRate(matrix,1);

  for (int n=1; n<int(matrix.getDimension()); n++) {
    double sqrtn=sqrt(n);
    dcomp offdiag(matrix(n)(n-1));
    averages(1)+=sqrtn*real(offdiag);
    averages(2)+=sqrtn*imag(offdiag);
  }

  return averages;

}


const Tridiagonal::Diagonal mainDiagonal(dcomp z, size_t cutoff);


basic::PumpedLossyModeIP::PumpedLossyModeIP(double delta, double kappa, dcomp eta, double nTh, size_t cutoff)
  : Free(cutoff,
         {
           RF{"kappa*(2*nTh+1)",kappa*(2*nTh+1),cutoff},
           RF{"delta",delta,1}
         },
         {CF{"eta",eta,sqrt(cutoff)}}),
    FreeExact<false>(cutoff),
    TridiagonalHamiltonian<1,true>(furnishWithFreqs(tridiagPlusHC_overI(conj(eta)*aop(cutoff)),
                                                    mainDiagonal(dcomp(kappa*(2*nTh+1),-delta),cutoff))),
    ElementLiouvilleanStrategies<1,2>(JumpStrategies(bind(aJump   ,_1,kappa*(nTh+1)),
                                                     bind(aDagJump,_1,kappa* nTh   )),
                                      JumpRateStrategies(bind(aJumpRate   ,_1,kappa*(nTh+1)),
                                                         bind(aDagJumpRate,_1,kappa* nTh   )),
                                      "Mode",{"photon loss","photon absorption"}),
    ElementAveraged<1>("PumpedLossyMode",{"<number operator>","real(<ladder operator>)","imag(\")"}),
    z_(kappa*(2*nTh+1),-delta)
{}


void basic::PumpedLossyModeIP::updateU(OneTime dtDid) const
{
  getDiagonal()=exp(-z_*(dtDid*blitz::tensor::i));
}


// PumpedLossyModeIP::average_v exactly the same as PumpedLossyMode::average_v above
auto basic::PumpedLossyModeIP::average_v(NoTime, const LazyDensityOperator& matrix) const -> const Averages
{
  auto averages(initializedAverages());

  averages(0)=aJumpRate(matrix,1);

  for (int n=1; n<int(matrix.getDimension()); n++) {
    double sqrtn=sqrt(n);
    dcomp offdiag(matrix(n)(n-1));
    averages(1)+=sqrtn*real(offdiag);
    averages(2)+=sqrtn*imag(offdiag);
  }

  return averages;

}


hierarchical::ModeBase::ModeBase(double kappa, double nTh, size_t cutoff)
  : Free(cutoff,{RF{"kappa*(2*nTh+1)",kappa*(2*nTh+1),cutoff}}),
    ElementLiouvillean<1,2>("Mode",{"photon loss","photon absorption"}),
    ElementAveraged<1>("PumpedLossyMode",{"<number operator>","real(<ladder operator>)","imag(\")"}),
    kappa_(kappa), nTh_(nTh)
{}


void hierarchical::ModeBase::doActWithJ(NoTime, StateVectorLow& psi, LindbladNo<0>) const
{
  aJump(psi,kappa_*(nTh_+1));
}

void hierarchical::ModeBase::doActWithJ(NoTime, StateVectorLow& psi, LindbladNo<1>) const
{
  aDagJump(psi,kappa_*nTh_);
}

double hierarchical::ModeBase::rate(NoTime, const LazyDensityOperator& matrix, LindbladNo<0>) const
{
  return aJumpRate(matrix,kappa_*(nTh_+1));
}

double hierarchical::ModeBase::rate(NoTime, const LazyDensityOperator& matrix, LindbladNo<1>) const
{
  return aDagJumpRate(matrix,kappa_*nTh_);
}

// ModeBase::average_v exactly the same as PumpedLossyMode::average_v above
auto hierarchical::ModeBase::average_v (NoTime, const hierarchical::ModeBase::LazyDensityOperator& matrix) const -> const Averages
{
  auto averages(initializedAverages());

  averages(0)=aJumpRate(matrix,1);

  for (int n=1; n<int(matrix.getDimension()); n++) {
    double sqrtn=sqrt(n);
    dcomp offdiag(matrix(n)(n-1));
    averages(1)+=sqrtn*real(offdiag);
    averages(2)+=sqrtn*imag(offdiag);
  }

  return averages;

}


hierarchical::PumpedLossyMode::PumpedLossyMode(double delta, double kappa, dcomp eta, double nTh, size_t cutoff)
  : ModeBase(kappa,nTh,cutoff),
    TridiagonalHamiltonian<1,false>(dcomp(-kappa*(2*nTh+1),delta)*nop(cutoff)
                                    +
                                    tridiagPlusHC_overI(conj(eta)*aop(cutoff)))
{}


hierarchical::PumpedLossyModeIP::PumpedLossyModeIP(double delta, double kappa, dcomp eta, double nTh, size_t cutoff)
  : ModeBase(kappa,nTh,cutoff),
    FreeExact<false>(cutoff),
    TridiagonalHamiltonian<1,true>(furnishWithFreqs(tridiagPlusHC_overI(conj(eta)*aop(cutoff)),
                                                    mainDiagonal(dcomp(kappa*(2*nTh+1),-delta),cutoff))),
    z_(kappa*(2*nTh+1),-delta)
{}


// PumpedLossyModeIP::updateU exactly the same as above
void hierarchical::PumpedLossyModeIP::updateU(OneTime dtDid) const
{
  getDiagonal()=exp(-z_*(dtDid*blitz::tensor::i));
}
