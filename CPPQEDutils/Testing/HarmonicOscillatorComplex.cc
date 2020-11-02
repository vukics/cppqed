// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ContainerIO.h"
#include "DrivenDampedHarmonicOscillator.h"
#include "Simulated.h"


typedef CArray<1> Array;

/*
  y(0) y
  y(1) dy/dt
*/


int main(int argc, char* argv[])
{
  ParameterTable p;

  Pars pt(p);

  double & omega=p.add("O","Driving frequency",5.),
         & gamma=p.add("G","Damping rate"     ,0.1);

  dcomp &    yinit=p.add(   "yinit"," y   initial condition",dcomp(10,5)),
        & dydtinit=p.add("dydtinit","dydt initial condition",dcomp(-1,1));

  // Parameter finalization
  update(p,argc,argv,"--");

  pt.T=70.; pt.dc=0; pt.Dt=0.1;

  // Initial condition
  Array y(2); y=yinit,dydtinit;

  Simulated<Array> S(y,
                     [=](double tau, const Array& y, Array& dydt)
                     {
                       dydt(0)=y(1);
                       dydt(1)=exp(DCOMP_I*omega*tau)-2*gamma*y(1)-y(0);
                     },
                     .1/std::max(1.,std::max(omega,gamma)),
                     pt);

  auto streamedArray=run(S,pt,false,true);

  ddho::Ptr oscillator(ddho::make(gamma,omega,yinit,dydtinit,0));

  // std::cout<<streamedArray;
  auto size=int(streamedArray.size());
  
  auto averageDt=0., ampDev=0., ampDerivDev=0.;
  
  for ( auto s=streamedArray.begin(); s!=streamedArray.end(); averageDt+=std::get<1>(*s++) ) {
    double time=std::get<0>(*s);
    
    auto amp=std::get<2>(*s)(0),
         ampDeriv=std::get<2>(*s)(1),
         exactAmp=oscillator->amp(time),
         exactAmpDeriv=oscillator->ampDeriv(time);    

    ampDev+=mathutils::relativeDeviation(amp,exactAmp);
    ampDerivDev+=mathutils::relativeDeviation(ampDeriv,exactAmpDeriv);

    // std::cout<<time<<"\t"<<std::get<1>(s)<<"\t"<<std::get<2>(s)(0)<<"\t"<<oscillator->amp(time)<<"\t"<<std::get<2>(s)(1)<<"\t"<<oscillator->ampDeriv(time)<<std::endl;
  }

  // std::cout<<averageDt/size<<"\t"<<ampDev/size<<"\t"<<ampDerivDev/size<<std::endl;
  
  return !(ampDev/size<2e-7 && ampDerivDev/size<2e-7);
  
}
