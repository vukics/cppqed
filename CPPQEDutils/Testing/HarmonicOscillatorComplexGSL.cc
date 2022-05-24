// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "DrivenDampedHarmonicOscillator.h"
#include "ODE_GSL.h"
#include "Simulated.h"


typedef CArray<1> Array;

/*
  y(0) y
  y(1) dy/dt
*/


int main(int, char**)
{
  const double omega{5.}, gamma{0.1};

  const dcomp yinit{10,5}, dydtinit{-1,1};

  // Initial condition
  Array y(2); y=yinit,dydtinit;

  auto S{simulated::make<ODE_EngineGSL<Array>>(y,
    [=](const Array& y, Array& dydt, double tau)
    {
      dydt(0)=y(1);
      dydt(1)=exp(DCOMP_I*omega*tau)-2*gamma*y(1)-y(0);
    },
    {"coordinate","velocity"},
    .1/std::max(1.,std::max(omega,gamma)),
    0,1e-6,1e-18)};

  auto streamedArray=run(S,70.,1,0,"","",6,false,false,"",false,true,trajectory::autostopHandlerNoOp);
  
  auto oscillator{ddho::make(gamma,omega,yinit,dydtinit,0)};

  auto size=int(streamedArray.size());
  
  auto averageDt=0., ampDev=0., ampDerivDev=0.;
  
  for ( auto s=streamedArray.begin(); s!=streamedArray.end(); averageDt+=std::get<1>(*s++) ) {
    double time=std::get<0>(*s);
    
    auto amp=std::get<2>(*s)(0),
         ampDeriv=std::get<2>(*s)(1),
         exactAmp=oscillator->amp(time),
         exactAmpDeriv=oscillator->ampDeriv(time);

    ampDev+=relativeDeviation(amp,exactAmp);
    ampDerivDev+=relativeDeviation(ampDeriv,exactAmpDeriv);

    // std::cout<<time<<"\t"<<std::get<1>(s)<<"\t"<<std::get<2>(s)(0)<<"\t"<<oscillator->amp(time)<<"\t"<<std::get<2>(s)(1)<<"\t"<<oscillator->ampDeriv(time)<<std::endl;
  }

  std::cout<<averageDt/size<<"\t"<<ampDev/size<<"\t"<<ampDerivDev/size<<std::endl;
  
  return !(ampDev/size<2e-6 && ampDerivDev/size<8e-7);
  
}
