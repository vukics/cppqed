#include "Evolution.h"
#include "Mode.h"
#include "MultiLevel.h"

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

  MultiLevelBase<NL>::Ptr plml(makePumpedLossyMultiLevelSch(pml,multilevel::DiagonalDO("Lambda atom",NL)));

  mode::Ptr mode(mode::make(pplm,QMP_IP));

  MLJC<NL,Couplings> mljc(plml,mode,pmljc);

  evolve<tmptools::Vector<0> >(psi,
			       binary::make(mljc),
			       pe);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}

}

