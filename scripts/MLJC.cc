#include "EvolutionBinary.h"
#include "Mode.h"
#include "MultiLevel.h"

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

  evolve<tmptools::Vector<0> >
    (psi,
     binary::make(MLJC<NL,Couplings>(makePumpedLossyMultiLevelSch(pml,multilevel::ReducedDensityOperator("Atom",NL)),
				     mode::make(pplm,QMP_IP),pmljc)),
     pe);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}

}