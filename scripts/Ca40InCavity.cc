#include "EvolutionHigh.h"

#include "MultiLevel.h"
#include "ParsMode.h"
#include "ParsMultiLevel.h"

#include "MLJC.h"
#include "ParsMLJC.h"

#include "BinarySystem.h"

#include "StateVector.h"

#include <boost/fusion/sequence/intrinsic/at.hpp>


using namespace std;

const int NL=8;

using multilevel::result_of::make_vector;
using multilevel::Pump; using multilevel::Decay;

typedef multilevel::RealLevelsMF<NL>::type Levels;

typedef make_vector<Pump<1,2>,Pump<0,3>/*,Pump<4,2>,Pump<6,2>,Pump<5,3>,Pump<7,3>*/ >::type Pumps;
typedef make_vector<Decay<0,2>,Decay<1,2>,Decay<0,3>,Decay<1,3>,Decay<4,2>,Decay<5,2>,Decay<6,2>,Decay<5,3>,Decay<6,3>,Decay<7,3> >::type Decays;

typedef make_vector<Coupling<4,2>,Coupling<5,2>,Coupling<6,2>,Coupling<5,3>,Coupling<6,3>,Coupling<7,3> >::type Couplings;


typedef quantumdata::StateVector<2> StateVector;

using boost::fusion::at_c;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  multilevel::ParsPumpedLossy<NL,Pumps,Decays> pml (p);
  mode::      ParsPumpedLossy                  pplm(p); 

  mljc::Pars<Couplings> pmljc(p);

  at_c<0>(pml.etas).set(30);
  at_c<1>(pml.etas).set(30);

  double gammaPS=10.48, gammaPD=.73, sqrt13=sqrt(1./3), sqrt23=sqrt(2./3), sqrt12=sqrt(1./2), sqrt16=sqrt(1./6);

  at_c<0>(pml.gammas).set(sqrt13*gammaPS);
  at_c<1>(pml.gammas).set(sqrt23*gammaPS);
  at_c<2>(pml.gammas).set(sqrt23*gammaPS);
  at_c<3>(pml.gammas).set(sqrt13*gammaPS);

  at_c<4>(pml.gammas).set(sqrt12*gammaPD);
  at_c<5>(pml.gammas).set(sqrt13*gammaPD);
  at_c<6>(pml.gammas).set(sqrt16*gammaPD);
  at_c<7>(pml.gammas).set(sqrt16*gammaPD);
  at_c<8>(pml.gammas).set(sqrt13*gammaPD);
  at_c<9>(pml.gammas).set(sqrt12*gammaPD);

  double g=1.4;

  at_c<0>(pmljc.gs).set(sqrt12*g);
  at_c<1>(pmljc.gs).set(sqrt13*g);
  at_c<2>(pmljc.gs).set(sqrt16*g);
  at_c<3>(pmljc.gs).set(sqrt16*g);
  at_c<4>(pmljc.gs).set(sqrt13*g);
  at_c<5>(pmljc.gs).set(sqrt12*g);

  double deltaEta=-318.803, deltaG=-315.12;

  pml.deltas(0)=0;
  
  pplm.kappa=.055;
  pplm.delta=deltaG-deltaEta;


  ParsEvolution pe(p); // Driver Parameters

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  multilevel::StateVector psiML(NL); psiML()=0; psiML()(1)=1;

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
