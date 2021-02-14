// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Simulated.h"

using namespace std;

typedef blitz::Array<double,1> DA1R;

/*
  y(0) Re y
  y(1) Im y
  y(2) Re dy/dt
  y(3) Im dy/dt
*/


int main(int argc, char* argv[])
{
  ParameterTable p;

  Pars pt(p);

  double& omega=p.add("O","Driving frequency",1.);
  double& gamma=p.add("G","Damping rate"     ,1.);

  dcomp&    yinit=p.add(   "yinit"," y   initial condition",dcomp( 1,-1));
  dcomp& dydtinit=p.add("dydtinit","dydt initial condition",dcomp(-1, 1));

  // Parameter finalization
  update(p,argc,argv,"--");

  if (pt.T<0) pt.T=10./min(1.,min(omega,gamma));
  // Note: 1.0 is also an existing frequency in the system, which defines the unit of time 

  DA1R y(4); y=yinit.real(),yinit.imag(),dydtinit.real(),dydtinit.imag();

  auto S{simulated::make(y,
    [=](const DA1R& y, DA1R& dydt, double tau) {
      dydt(0)=y(2); dydt(1)=y(3); 
      dydt(2)=cos(omega*tau)-2*gamma*y(2)-y(0); dydt(3)=sin(omega*tau)-2*gamma*y(3)-y(1);
    },
    .1/max(1.,max(omega,gamma)), 
    {"Re{coordinate}","Im{coordinate}","Re{velocity}","Im{velocity}"},
    1e-6,1e-18)};
  
  run(S,pt);

}
