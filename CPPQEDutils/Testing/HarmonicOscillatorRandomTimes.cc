// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Simulated.h"

#include "DrivenDampedHarmonicOscillator.h"
#include "MathExtensions.h"
#include "Random.h"

typedef CArray<1> Array;

/*
  y(0) y
  y(1) dy/dt
*/


int main(int, char**)
{
  const double omega{5.}, gamma{0.1}, T{70.};

  const dcomp yinit{10,5}, dydtinit{-1,1};

  const size_t nSamples=1000;

  // Initial condition
  Array y(2); y=yinit,dydtinit;

  Simulated<Array> s(y,
                     [=](double tau, const Array& y, Array& dydt)
                     {
                       dydt(0)=y(1);
                       dydt(1)=exp(DCOMP_I*omega*tau)-2*gamma*y(1)-y(0);
                     },
                     .1/std::max(1.,std::max(omega,gamma)),
                     0,1e-6,1e-18);
  
  // s.streamParameters(std::cout);

  cpputils::GSL_RandomEngine re{1001};

  auto ampDev=0., ampDerivDev=0.;

  auto oscillator{ddho::make(gamma,omega,yinit,dydtinit,0)};
  
  for (size_t n=nSamples; n>0; n--) {
    auto time=std::uniform_real_distribution(0.,T)(re);
    evolved::evolveTo(s,time);
    
    ampDev+=mathutils::relativeDeviation(y(0),oscillator->amp(time));
    ampDerivDev+=mathutils::relativeDeviation(y(1),oscillator->ampDeriv(time));

  }

  // std::cout<<ampDev/nSamples<<"\t"<<ampDerivDev/nSamples<<std::endl;

  return !(ampDev/nSamples<1e-3 && ampDerivDev/nSamples<1e-3);
  
}
