// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"
#include "MLJC.h"


using namespace std;
using namespace multilevel;

const int NL=3; // NL stands for "Number of Levels"


using Drives=hana::tuple< Pump<0,2>, Pump <1,2> >;

using Decays=hana::tuple< Decay<0,2>, Decay<1,2> >;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  ParsDrivenDissipative<NL,Drives,Decays> pml(p);

  evolution::Pars<> pe(p); // Driver Parameters

  pml.etas  =hana::make_tuple(Drive<0,2>{30,40},Drive<1,2>{-60,40});
  
  pml.gammas=hana::make_tuple(Decay<0,2>{10},Decay<1,2>{5});

  pml.deltas=RealPerLevel<3>(0,0,-1e5);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  StateVector psi{NL}; psi(0)=1; psi(1)=1; psi.renorm();

  evolve(psi,
         multilevel::makeDrivenDissipativeSch(pml,"Lambda atom",true),
         pe);



}
