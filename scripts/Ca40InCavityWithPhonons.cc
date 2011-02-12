#include "EvolutionHigh.h"

#include "MultiLevel.h"
#include "ParsMode.h"
#include "ParsMultiLevel.h"

#include "MLJC.h"
#include "ParsMLJC.h"

#include "Composite.h"

#include <boost/fusion/sequence/intrinsic/at.hpp>


using namespace std;

const int NL=8;

using multilevel::result_of::make_vector;
using multilevel::Pump; using multilevel::Decay;

typedef multilevel::RealLevelsMF<NL>::type Levels;

typedef make_vector<Pump    <1,2>,Pump    <0,3>/*,Pump<4,2>,Pump<6,2>,Pump<5,3>,Pump<7,3>*/ >::type Pumps          ;
typedef make_vector<Coupling<1,2>,Coupling<0,3>                                             >::type PhononCouplings;

typedef make_vector<Decay<0,2>,Decay<1,2>,Decay<0,3>,Decay<1,3>,Decay<4,2>,Decay<5,2>,Decay<6,2>,Decay<5,3>,Decay<6,3>,Decay<7,3> >::type Decays;

typedef make_vector<Coupling<4,2>,Coupling<5,2>,Coupling<6,2>,Coupling<5,3>,Coupling<6,3>,Coupling<7,3> >::type Couplings;


typedef quantumdata::StateVector<3> StateVector;

using boost::fusion::at_c;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  // Clebsches:
  const double sqrt13=sqrt(1./3), sqrt23=sqrt(2./3), sqrt12=sqrt(1./2), sqrt16=sqrt(1./6);

  ParameterTable pDummy, pActual;

  multilevel::ParsPumpedLossy<NL,Pumps,Decays> pml  (pDummy);
  mode::      Pars                             pphon(pDummy,"P");
  mode::      ParsPumpedLossy                  pplm (pDummy,"C"); 

  mljc::Pars<PhononCouplings> pPhononCouplings(pDummy,"P");
  mljc::Pars<      Couplings> pmljc           (pDummy,"C");

  {
    // These can be hardcoded here:

    const double gammaPS=10.48, gammaPD=.73;

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
  }

  dcomp
    &eta =pActual.add("eta","Pumping laser Rabi frequency",dcomp(30,0)),
    &etaC=pActual.add("etaC","Cavity pump Rabi frequency",dcomp(0,0));

  double 
    &phononCouplingFactor=pActual.add("phononCouplingFactor","Phonon coupling strength: -i*eta*phononCouplingFactor",.12),
    &g                   =pActual.add("g","Atom-cavity coupling: Clebsches*g",1.4),
    &deltaEta            =pActual.add("deltaEta","Pump detuning",-318.803), 
    &deltaG              =pActual.add("deltaG","Cavity detuning",-315.12 ),
    &omegaPhonon         =pActual.add("omegaPhonon","Phonon frequency",1.);

  unsigned
    &phononCutoff=pActual.add("phononCutoff","Phonon mode cutoff",20u),
    &cavityCutoff=pActual.add("cavityCutoff","Cavity mode cutoff",3u);

  ParsEvolution pe(pActual); // Driver Parameters

  // Parameter finalization
  update(pActual,argc,argv,"--");

  {
    // Enforcing constraints:

    at_c<0>(pml.etas).set(eta);
    at_c<1>(pml.etas).set(eta);

    at_c<0>(pPhononCouplings.gs).set(-DCOMP_I*phononCouplingFactor*at_c<0>(pml.etas).get());
    at_c<1>(pPhononCouplings.gs).set(-DCOMP_I*phononCouplingFactor*at_c<1>(pml.etas).get());

    at_c<0>(pmljc.gs).set(sqrt12*g);
    at_c<1>(pmljc.gs).set(sqrt13*g);
    at_c<2>(pmljc.gs).set(sqrt16*g);
    at_c<3>(pmljc.gs).set(sqrt16*g);
    at_c<4>(pmljc.gs).set(sqrt13*g);
    at_c<5>(pmljc.gs).set(sqrt12*g);

    pml.deltas(0)=0;
  
    pplm.kappa=.055;
    pplm.delta=deltaG-deltaEta;

    pphon.cutoff=phononCutoff;
    pplm .cutoff=cavityCutoff;

    pphon.delta=-omegaPhonon;

    pplm.eta=etaC;
  }

  // ****** ****** ****** ****** ****** ******

  multilevel::StateVector psiML(NL); psiML()=0; psiML()(1)=1;

  StateVector psi(psiML*mode::init(pphon)*mode::init(pplm));

  // The free components

  PumpedLossyMultiLevelSch<NL,Pumps,Decays> atomInner(pml.deltas,pml.etas,pml.gammas);

  Mode<> phonon(pphon);
  
  mode::SmartPtr cavityMode(mode::maker(pplm,QMP_IP));

  // The interaction components


  MLJC<NL,PhononCouplings> ia01(atomInner,phonon    ,pPhononCouplings);

  MLJC<NL,      Couplings> ia02(atomInner,cavityMode,pmljc           );

  evolve(psi,
	 makeComposite(Act<0,1>(ia01),Act<0,2>(ia02)),
	 pe,
	 tmptools::Vector<0>());

  } catch (ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}

}

