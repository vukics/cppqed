// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Simulated.h"

using std::max, std::min;

typedef blitz::Array<dcomp,1> Array;

/*
  y(0) y
  y(1) dy/dt
*/


int main(int argc, char* argv[])
{
  ParameterTable p;

  Pars pt(p);

  double & omega=p.add("O","Driving frequency",1.),
         & gamma=p.add("G","Damping rate"     ,1.);

  dcomp &    yinit=p.add(   "yinit"," y   initial condition",dcomp(1,-1)/*,dcomp( 3.45263,-2.0746 )*/),
        & dydtinit=p.add("dydtinit","dydt initial condition",dcomp(-1,1)/*,dcomp(-4.83065, 1.16527)*/);

  // Parameter finalization
  update(p,argc,argv,"--");

  if (pt.T<0) pt.T=10./min(1.,min(omega,gamma));
  // Note: 1.0 is also an existing frequency in the system, which defines the unit of time 

  // Initial condition
  Array y(2); y=yinit,dydtinit;

  Simulated<Array> S(y,
                     [=](double tau, const Array& y, Array& dydt)
                     {
                       dydt(0)=y(1);
                       dydt(1)=exp(DCOMP_I*omega*tau)-2*gamma*y(1)-y(0);
                     },
                     .1/max(1.,max(omega,gamma)),
                     pt);

  run(S,pt);

}
