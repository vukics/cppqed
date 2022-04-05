// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "DrivenDampedHarmonicOscillator.h"

#include "MathExtensions.h"

#include <iostream>


dcomp cppqedutils::ddho::_::c(double t) const
{
  return exp(DCOMP_I*omega_*t)/(1.+2.*DCOMP_I*gamma_*omega_-sqr(omega_));
}


cppqedutils::ddho::_::_(double gamma, double omega, Matrix m, dcomp ampTI, dcomp ampDerivTI, double tInit)
  : gamma_(gamma), omega_(omega), a_(m.fullPivLu().solve(Vector(ampTI-c(tInit),ampDerivTI-DCOMP_I*omega*c(tInit))))
{}


namespace {


using namespace cppqedutils;


class DDHO : private std::pair<dcomp,dcomp>, public ddho::_
{
public:
  DDHO(double gamma, double omega, dcomp ampTI, dcomp ampDerivTI, double tInit=0)
    : std::pair<dcomp,dcomp>(-gamma+sqrt(dcomp(sqr(gamma)-1)),-gamma-sqrt(dcomp(sqr(gamma)-1))),
      _(gamma,omega,
        [=,this]() {ddho::Matrix m; m << exp(first*tInit), exp(second*tInit), first*exp(first*tInit), second*exp(second*tInit); return m;} (),
        ampTI,ampDerivTI,tInit),
      omega1_(first), omega2_(second) {}

  dcomp amp     (double t) const override {return         exp(omega1_*t)*a_(0)+        exp(omega2_*t)*a_(1)+               c(t);}
  dcomp ampDeriv(double t) const override {return omega1_*exp(omega1_*t)*a_(0)+omega2_*exp(omega2_*t)*a_(1)+DCOMP_I*omega_*c(t);}

private:
  const dcomp& omega1_, omega2_;

};


class DDHO_Critical : public ddho::_
{
public:
  DDHO_Critical(double omega, dcomp ampTI, dcomp ampDerivTI, double tInit=0)
    : _(1.,omega,
        [=]() {ddho::Matrix m; m << exp(-tInit), tInit*exp(-tInit), -exp(-tInit),-(tInit-1)*exp(-tInit); return m;} (),
        ampTI,ampDerivTI,tInit)
  {}

  dcomp amp     (double t) const override {return  exp(-t)*a_(0)+ t   *exp(-t)*a_(1)+               c(t);}
  dcomp ampDeriv(double t) const override {return -exp(-t)*a_(0)-(t-1)*exp(-t)*a_(1)+DCOMP_I*omega_*c(t);}

};


}


cppqedutils::ddho::Ptr cppqedutils::ddho::make(double gamma, double omega, dcomp ampTI, dcomp ampDerivTI, double tInit)
{
  if (gamma==1.)
    return std::make_shared<DDHO_Critical>(omega,ampTI,ampDerivTI,tInit);
  else
    return std::make_shared<DDHO>(gamma,omega,ampTI,ampDerivTI,tInit);
}

