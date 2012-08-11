#include "Evolution.h"
#include "Mode.h"
#include "MCWF.h"
#include "MultiLevel.h"

#include "BinarySystem.h"

#include "StateVector.h"

using namespace std;

const int NL=5;

using multilevel::result_of::make_vector;
using multilevel::Pump; using multilevel::Decay;

typedef multilevel::RealLevelsMF<NL>::type Levels;

typedef make_vector<Pump<3,2>,Pump<3,1> >::type Pumps;
typedef make_vector<Decay<3,2>,Decay<2,1> >::type Decays;

typedef make_vector<Coupling<4,2>,Coupling<3,1>,Coupling<4,1> >::type Couplings;


typedef quantumdata::StateVector<2> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  multilevel::ParsPumpedLossy<NL,Pumps,Decays> pml (p);
  mode::      ParsPumpedLossy                  pplm(p); 

  mljc::Pars<Couplings> pmljc(p);

  ParsEvolution pe(p); // Driver Parameters

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  multilevel::StateVector psiML(NL); psiML()(0)=1;

  StateVector psi(psiML*mode::init(pplm));

  evolve(psi,
	 binary::make(MLJC<NL,Couplings>(makePumpedLossyMultiLevelSch(pml,multilevel::DiagonalDO("Atom",NL)),mode::make(pplm,QMP_IP),pmljc)),
	 pe,
	 tmptools::Vector<0>());

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}

}


/*

Note that the following snippet yields undefined behaviour (clearly!)


1. MLJC<NL,Couplings> mljc(*makePumpedLossyMultiLevelSch(pml,multilevel::DiagonalDO(NL)),*mode::make(pplm,QMP_IP),gs);

2. evolve(psi,BinarySystem(mljc),pe);


In Line 1 TEMPORARY SmartPtrs are created with the maker functions and
are COPIED into mljc::Interaction as plain ptrs. "Between" Lines 1 and
2 the SmartPtrs are destructed and the ptrs deleted. So, they become
dangling ptrs in mljc!


*/
