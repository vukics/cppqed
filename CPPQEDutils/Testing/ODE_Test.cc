// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/*
  
  A rather artificial example based on the observation that

  exp(exp(Zt))

  is a particular solution of the equation

  D^2(y) - D(y)^2 - Z D(y) + Z^2 exp(2Zt) ( y^2 - y ) = 0

  for initial condition

  y(0) = e
  p(0) = Z e

 */

#include "ComplexExtensions.h"
#include "Simulated.h"

#include <numbers>

using namespace std;
using namespace trajectory;

using Array=std::array<dcomp,2>;


int main(int , char**)
{
  dcomp Z(-1.,10);

  Array y{numbers::e,Z*numbers::e};

  auto S(simulated::make<ODE_EngineBoost<Array>>(y,[=](const Array& yA, Array& dydtA, double tau) {
    dydtA[0]=yA[1];
    dydtA[1]=sqr(yA[1])+Z*yA[1]+sqr(Z)*exp(2.*Z*tau)*(yA[0]-sqr(yA[0]));
  },{"Re{y}","Im{y}"},trajectory::initialTimeStep(abs(Z)),
  0, //logLevel
  1e-6,1e-18));

  auto streamedArray=run<RunLengthType::T_MODE,StreamFreqType::DT_MODE>(S,5.,0.01,0,"","",6,false,false,false,true,autostopHandlerNoOp);

  double yDev=0;
  
  for ( auto s : streamedArray ) yDev+=relativeDeviation(get<2>(s)[0],exp(exp(Z*get<0>(s))));
  
  // for ( auto s : streamedArray ) cout<<get<0>(s)<<"\t"<<get<2>(s)(0).real()<<"\t"<<get<2>(s)(0).imag()<<endl;
  
  return !(yDev/streamedArray.size()<1.5e-8);
 
}
