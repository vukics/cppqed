// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "DrivenDampedHarmonicOscillator.h"

#include "Simulated.h"

using namespace std       ;
using namespace trajectory;
using namespace mathutils ;
using namespace parameters;


int main(int argc, char* argv[])
{
  ParameterTable p;

  Pars pt(p);

  double
    &omega=p.add("O","Driving frequency",1.),
    &gamma=p.add("G","Damping rate"     ,1.);
  
  // Parameter finalization
  update(p,argc,argv,"--");

  if (pt.T<0) pt.T=10./min(1.,min(omega,gamma));

  double dtInit=.1/max(1.,max(omega,gamma));  

  ddho::Ptr oscillator(ddho::make(gamma,omega,dcomp(1.,-1.),dcomp(-1.,1.),0));

  for (double t=0.; t<pt.T; t+=dtInit)
    cout<<t<<' '<<dtInit<<' '<<oscillator->amp(t)<<' '<<oscillator->ampDeriv(t)<<endl;

}
