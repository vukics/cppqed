// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Simulated.h"

#include <boost/bind.hpp>

using std::max;

typedef blitz::Array<dcomp,1> Array;

/*
  y(0) y
  y(1) dy/dt
*/

void derivs(double tau, const Array& y, Array& dydt, double omega, double gamma)
{
  dydt(0)=y(1);
  dydt(1)=exp(DCOMP_I*omega*tau)-2*gamma*y(1)-y(0);
}


int main(int argc, char* argv[])
{
  ParameterTable pt;

  Pars p(pt);

  double & omega=pt.add("O","Driving frequency",1.),
         & gamma=pt.add("G","Damping rate"     ,1.);

  dcomp &    yinit=pt.add(   "yinit"," y   initial condition",dcomp( 3.45263,-2.0746 )),
        & dydtinit=pt.add("dydtinit","dydt initial condition",dcomp(-4.83065, 1.16527));

  // Parameter finalization
  update(pt,argc,argv,"--");

  // Initial condition
  Array y(2); y=yinit,dydtinit;

  Simulated<Array> S(y,
                     bind(derivs,_1,_2,_3,omega,gamma),
                     .1/max(1.,max(omega,gamma)),
                     p);

  run(S,p);

}
