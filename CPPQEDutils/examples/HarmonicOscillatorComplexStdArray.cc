// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Simulated.h"

using namespace std;

using Array = std::array<dcomp,2>;

/*
  y(0) y
  y(1) dy/dt
*/


int main(int argc, char* argv[])
{
  auto op{optionParser()};

  simulated::Pars<> pt(op);

  double omega, gamma;
  dcomp yinit, dydtinit;
  
  add(add(add(add(op,
      "dydtinit","dydt initial condition",dcomp(-1,1)/*,dcomp(-4.83065, 1.16527)*/,&dydtinit),
      "yinit","y initial condition",dcomp(1,-1)/*,dcomp( 3.45263,-2.0746 )*/,&yinit),
      "G","Damping rate",1.,&gamma),
      "O","Driving frequency",1.,&omega);

  parse(op,argc, argv);
  
  if (pt.T<0) pt.T=10./min(1.,min(omega,gamma));
  // Note: 1.0 is also an existing frequency in the system, which defines the unit of time 

  // Initial condition
  Array y{yinit,dydtinit};

  auto S{simulated::makeBoost(y,
    [=](const Array& y, Array& dydt, double tau)
    {
      dydt[0]=y[1];
      dydt[1]=exp(1i*omega*tau)-2*gamma*y[1]-y[0];
    },
    {"complex coordinate","complex velocity"},
    {{"G",gamma},{"O",omega}},
    .1/max(1.,max(omega,gamma)),pt)};

  run(S,pt,trajectory::observerNoOp);
/*
  auto streamedArray=run(S,pt);
  
  for (const auto & t : streamedArray) std::cout<<std::get<0>(t)<<"\t"<<std::get<1>(t)<<"\t"<<std::get<2>(t);
*/
}
