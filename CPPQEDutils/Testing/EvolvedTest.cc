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

using namespace std;
using namespace mathutils;
using namespace trajectory;

typedef CArray<1> Array;


int main(int , char**)
{
  dcomp Z(-1.,10);

  Array y(2);

  y(0)=  EULER;
  y(1)=Z*EULER;

  Simulated<Array> S(y,[=](double tau, const Array& yA, Array& dydtA) {
                        dcomp y{yA(0)}, p{yA(1)};
                        dydtA(0)=p;
                        dydtA(1)=sqr(p)+Z*p+sqr(Z)*exp(2.*Z*tau)*(y-sqr(y));
                      },
                    trajectory::initialTimeStep(abs(Z)),0,1e-6,1e-18);

  auto streamedArray=run(static_cast<Trajectory<Array>&>(S),5.,0.01,0,string{},string{},6,false,false,string{},false,true,
                         AutostopHandlerNoOp<Array>{});

  double yDev=0;
  
  for ( auto s : streamedArray ) yDev+=relativeDeviation(std::get<2>(s)(0),exp(exp(Z*std::get<0>(s))));
  
  // for ( auto s : streamedArray ) std::cout<<std::get<0>(s)<<"\t"<<std::get<2>(s)(0).real()<<"\t"<<std::get<2>(s)(0).imag()<<std::endl;
  
  return !(yDev/streamedArray.size()<3e-8);
 
}
