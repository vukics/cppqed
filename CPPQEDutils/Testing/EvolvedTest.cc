// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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

#include <boost/bind.hpp>


using namespace std;
// using namespace blitzplusplus;
using namespace mathutils;
using namespace parameters;
using namespace trajectory;
/*
  y[0] Re y
  y[1] Im y
  y[2] Re dy/dt
  y[3] Im dy/dt
*/

// #include<cmath>
// #include<iostream>

typedef CArray<1> CA1D;
typedef DArray<1> DA1D;

void derivs(double tau, const CA1D& yA, CA1D& dydtA, const dcomp& Z)
{
  // static int count=0;

  dcomp
    y(yA(0)),
    p(yA(1));

  dydtA(0)=p;

  dydtA(1)=sqr(p)+Z*p+sqr(Z)*exp(2.*Z*tau)*(y-sqr(y));

  // cout<<count++<<endl;

}


int main(int , char**)
{
  ParameterTable p;

  Pars pt(p);

  dcomp Z(-1.,10);

  pt.T=5; pt.dc=0; pt.Dt=0.01;
  
  CA1D y(2);

  y(0)=  EULER;
  y(1)=Z*EULER;

  Simulated<CA1D> S(y,bind(derivs,_1,_2,_3,Z),trajectory::initialTimeStep(abs(Z)),pt);

  auto streamedArray=run(S,pt,false,true);

  double yDev=0;
  
  for ( auto s : streamedArray ) yDev+=relativeDeviation(std::get<2>(s)(0),exp(exp(Z*std::get<0>(s))));
  
  // for ( auto s : streamedArray ) std::cout<<std::get<0>(s)<<"\t"<<std::get<2>(s)(0).real()<<"\t"<<std::get<2>(s)(0).imag()<<std::endl;
  
  return !(yDev/streamedArray.size()<3e-8);
 
}
