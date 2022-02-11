// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"
#include "MLJC.h"


using namespace std;
using namespace multilevel;

const int NL=3; // NL stands for "Number of Levels"

typedef multilevel::result_of::make_vector<Pump <0,2>,Pump <1,2> >::type Pumps;

typedef multilevel::result_of::make_vector<Decay<0,2>,Decay<1,2> >::type Decays;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  ParsPumpedLossy<NL,Pumps,Decays> pml(p);

  evolution::Pars<> pe(p); // Driver Parameters

  pml.etas=make_vector(dcomp(30,40),dcomp(-60,40));
  pml.gammas=make_vector(10,5);

  pml.deltas=RealPerLevel<3>(0,0,-1e5);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  StateVector psi{NL}; psi(0)=1; psi(1)=1; psi.renorm();

  evolve(psi,
         multilevel::makePumpedLossySch(pml,"Lambda atom",true),
         pe);



}
