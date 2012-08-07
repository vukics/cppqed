#include "EvolutionHigh.h"

#include "MultiLevel.h"
#include "ParsMultiLevel.h"

#include "StateVector.h"

using namespace std;
using namespace multilevel;

const int NL=3; // NL stands for "Number of Levels"

typedef RealLevelsMF<NL>::type Levels;

typedef result_of::make_vector<Pump <0,2>,Pump <1,2> >::type Pumps;

typedef result_of::make_vector<Decay<0,2>,Decay<1,2> >::type Decays;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsPumpedLossy<NL,Pumps,Decays> pml(p);

  ParsEvolution pe(p); // Driver Parameters

  pml.etas=make_vector(dcomp(30,40),dcomp(-60,40));
  pml.gammas=make_vector(10,5);

  pml.deltas=Levels(0,0,-1e5);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  StateVector psi(NL); psi()(0)=1; psi()(1)=1; psi.renorm();

  evolve(psi,
	 makePumpedLossyMultiLevelSch(pml,DiagonalDO("Lambda atom",NL)),
	 pe);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}

}
