// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/*
  
  A rather artificial example based on the observation that

  exp(exp(Zt))

  is a particular solution of the equation

  D^2(y) - D(y)^2 - Z D(y) + Z^2 exp(2Zt) ( y^2 - y ) = 0

  for initial condition

  y(0) = e
  p(0) = Z e

 */

#include "Simulated.h"

#include "MathExtensions.h"

using namespace std;
using namespace trajectory;

typedef CArray<1> Array;


int main(int , char**)
{
  dcomp Z(-1.,10);

  Array y(2);

  y(0)=  EULER;
  y(1)=Z*EULER;

  auto S(simulated::makeBoost(y,[=](const Array& yA, Array& dydtA, double tau) {
    dcomp y{yA(0)}, p{yA(1)};
    dydtA(0)=p;
    dydtA(1)=sqr(p)+Z*p+sqr(Z)*exp(2.*Z*tau)*(y-sqr(y));
  },{"Re{y}","Im{y}"},trajectory::initialTimeStep(abs(Z)),
                              0, //logLevel
                              1e-6,1e-18));

  auto streamedArray=run(S,5.,0.01,0,string{},string{},6,false,false,string{},false,true,autostopHandlerNoOp);

  double yDev=0;
  
  for ( auto s : streamedArray ) yDev+=relativeDeviation(get<2>(s)(0),exp(exp(Z*get<0>(s))));
  
  // for ( auto s : streamedArray ) cout<<get<0>(s)<<"\t"<<get<2>(s)(0).real()<<"\t"<<get<2>(s)(0).imag()<<endl;
  
  return !(yDev/streamedArray.size()<2e-8);
 
}
