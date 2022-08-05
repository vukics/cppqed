// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "JaynesCummings.h"

#include "BinarySystem.h"


using namespace std;

typedef quantumdata::StateVector<2> StateVector;
typedef StateVector::Dimensions     Dimensions ;

int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  parameters::Table p;

  evolution::Pars<> pe(p); // Driver Parameters
  ParsPumpedLossyQbit pplqb(p); 
  ParsPumpedLossyMode pplm (p); 
  ParsJaynesCummings  pjc  (p); 

  pplqb.eta=dcomp(-.2, .3);
  pplm .eta=dcomp( .5,-.3);
  pjc  .g  =dcomp(  1,-.5);

  pplqb.gamma=20;
  pplqb.delta= 5;
  pplm .kappa=12;

  // Parameter finalization
  parameters::update(p,argc,argv,"--");

  QM_Picture qmp=QMP_IP;
  
  // ****** ****** ****** ****** ****** ******

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  mode::Ptr mode(mode::make(pplm ,qmp));

  JaynesCummings<true> jc(*qbit,*mode,pjc);

  BinarySystem system(jc);

  cout<<dynamic_cast<mode::Hamiltonian<true>*>(mode.get())->getFreqss().front()<<endl;
  cout<<dynamic_cast<qbit::Hamiltonian<true>*>(qbit.get())->getFreqss().front()<<endl;

  cout<<jc.getH_OverIs().front()<<endl;

  cout<<jc.getFreqss  ().front()<<endl;

  StateVector 
    psi   (qbit::state0()*mode::fock(3,pplm.cutoff)+qbit::state1()*mode::fock(8,pplm.cutoff)),
    dpsidt(psi.getDimensions());

  cout<<psi()<<endl;
  structure::Hamiltonian<2>::addContribution(&system,0,psi(),dpsidt(),0);
  cout<<dpsidt()<<endl;

  cout<<structure::Liouvillian<2>::probabilities(&system,psi)<<endl;


  structure::Liouvillian<2>::actWithJ(&system,psi(),0);
  cout<<psi()<<endl;

  structure::Liouvillian<2>::actWithJ(&system,psi(),1);
  cout<<psi()<<endl;


  /*()+qbit::state1()())*
  psi()(0,0)=1; psi()(1,0)=1;
  */
  // psi.renorm();

  evolve(psi,system,pe);

  } catch (const parameters::NamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
