#include "cppqedconfig.h"

#ifndef   DO_NOT_USE_FLENS


#include "DrivenDampedHarmonicOscillator.h"

#include "MathExtensions.h"

#include <boost/make_shared.hpp>



const dcomp DrivenDampedHarmonicOscillator::c(double t) const
{
  return exp(DCOMP_I*omega_*t)/(1.+2.*DCOMP_I*gamma_*omega_-mathutils::sqr(omega_));
}



const DrivenDampedHarmonicOscillator::Vector DrivenDampedHarmonicOscillator::calculateA(Matrix m, dcomp ampTI, dcomp ampDerivTI, dcomp cTInit, double omega)
{
  Vector rhs(2); rhs=ampTI-cTInit,ampDerivTI-DCOMP_I*omega*cTInit;

  {
    using namespace flens;
    DenseVector<Array<int> > pivots(2);
    trf(m,pivots);
    tri(m,pivots);
  }

  return m*rhs;

}


const DrivenDampedHarmonicOscillator::Matrix DrivenDampedHarmonicOscillator::makeMatrix(const dcomp& a0, const dcomp& a1, const dcomp& a2, const dcomp& a3)
{
  Matrix res(2,2);
  res=
    a0, a1,
    a2, a3;
  return res;
}


DrivenDampedHarmonicOscillator::DrivenDampedHarmonicOscillator(double gamma, double omega, const Matrix& m, dcomp ampTI, dcomp ampDerivTI, double tInit)
  : gamma_(gamma), omega_(omega), a_(calculateA(m,ampTI,ampDerivTI,c(tInit),omega))
{}





class DDHO : private std::pair<dcomp,dcomp>, public DrivenDampedHarmonicOscillator
{
public:
  DDHO(double gamma, double omega, dcomp ampTI, dcomp ampDerivTI, double tInit=0)
    : std::pair<dcomp,dcomp>(-gamma+sqrt(dcomp(mathutils::sqr(gamma)-1)),-gamma-sqrt(dcomp(mathutils::sqr(gamma)-1))),
      DrivenDampedHarmonicOscillator(gamma,omega,
				     makeMatrix(      exp(first*tInit),       exp(second*tInit),
						first*exp(first*tInit),second*exp(second*tInit) ),
				     ampTI,ampDerivTI,tInit),
      omega1_(first), omega2_(second) {}

  const dcomp amp     (double t) const {return         exp(omega1_*t)*a_(1)+        exp(omega2_*t)*a_(2)+               c(t);}
  const dcomp ampDeriv(double t) const {return omega1_*exp(omega1_*t)*a_(1)+omega2_*exp(omega2_*t)*a_(2)+DCOMP_I*omega_*c(t);}

private:
  const dcomp& omega1_, omega2_;

};


class DDHO_Critical : public DrivenDampedHarmonicOscillator
{
public:
  DDHO_Critical(double omega, dcomp ampTI, dcomp ampDerivTI, double tInit=0)
    : DrivenDampedHarmonicOscillator(1.,omega,
				     makeMatrix( exp(-tInit),  tInit   *exp(-tInit),
						-exp(-tInit),-(tInit-1)*exp(-tInit)),
				     ampTI,ampDerivTI,tInit)
  {}

  const dcomp amp     (double t) const {return  exp(-t)*a_(1)+ t   *exp(-t)*a_(2)+               c(t);}
  const dcomp ampDeriv(double t) const {return -exp(-t)*a_(1)-(t-1)*exp(-t)*a_(2)+DCOMP_I*omega_*c(t);}

};



const DDHO_Ptr makeDDHO(double gamma, double omega, dcomp ampTI, dcomp ampDerivTI, double tInit)
{
  if (gamma==1.)
    return boost::make_shared<DDHO_Critical>(omega,ampTI,ampDerivTI,tInit);
  else
    return boost::make_shared<DDHO>(gamma,omega,ampTI,ampDerivTI,tInit);
}


#endif // DO_NOT_USE_FLENS
