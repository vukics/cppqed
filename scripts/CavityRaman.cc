#include "EvolutionHigh.h"

#include "MultiLevel.h"
#include "ParsMode.h"
#include "ParsMultiLevel.h"

#include "MLJC.h"
#include "ParsMLJC.h"

#include "BinarySystem.h"

#include "StateVector.h"

using namespace std;

const int NL=3;

using multilevel::result_of::make_vector;
using multilevel::Pump; using multilevel::Decay;

typedef multilevel::RealLevelsMF<NL>::type Levels;

typedef make_vector<Pump <0,2>            >::type Pumps;
typedef make_vector<Decay<0,2>,Decay<1,2> >::type Decays;

typedef make_vector<Coupling<1,2> >::type Couplings;


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

  multilevel::StateVector psiML(NL); psiML()(0)=1; psiML()(1)=1; 

  StateVector psi(psiML*mode::init(pplm)); psi.renorm();

  multilevel::DiagonalDO averaged(NL);

  MultiLevelBase<NL>::SmartPtr plml(makePumpedLossyMultiLevelSch(pml,multilevel::DiagonalDO(NL)));

  mode::SmartPtr mode(mode::maker(pplm,QMP_IP));

  MLJC<NL,Couplings> mljc(plml,mode,pmljc);

  evolve(psi,
	 BinarySystem(mljc),
	 pe,
	 tmptools::Vector<0>());

  } catch (ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}

}

