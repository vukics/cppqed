// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// #include "Evolution.h"

#include "Qbit.h"

#include "Master.h"

#include <iostream>

using namespace qbit;
using namespace std;

void f(hamiltonian<1> auto) {}

int main(int argc, char* argv[])
{
  auto op{optionParser()};

  cppqedutils::trajectory::Pars<quantumtrajectory::master::Pars> pt(op);

  parse(op,argc, argv);

  dcomp z{.1,1.}, eta{-.1,.2};

  QuantumSystemDynamics qbit {
    std::array{structure::SystemFrequencyDescriptor{"z",z,1.},structure::SystemFrequencyDescriptor{"η",eta,1.}},
    std::array{loss(1.),gain(.1)},
    makeHamiltonianCollection<1>(diagonalH(z),offDiagonalH(eta)),
    expectationValues }; //= make({1,1},{1,1},{1,1},1,1,1);
  
  auto oe=::cppqedutils::ODE_EngineBoost<::quantumdata::DensityOperator<1>::StorageType>(1./(10.*highestFrequency(qbit.freqs)),1e-12,1e-30);

  quantumtrajectory::Master<1,decltype(qbit),decltype(oe)> m{qbit,quantumdata::DensityOperator<1>{{2}},oe};

  cppqedutils::run(m,pt,cppqedutils::trajectory::observerNoOp);

/*  for (size_t i=0; i<1e2; ++i) {
    step(m,1.); cppqedutils::trajectory::dataStreamerDefault(m,std::cerr<<getTime(m)<<" "<<getDtDid(m)<<"\t");
  }

  auto streamedArray=run<cppqedutils::trajectory::RunLengthType::T_MODE,cppqedutils::trajectory::StreamFreqType::DC_MODE>(
    m,70.,1,0,"","",6,cppqedutils::trajectory::StreamSwitch{"1111"},false,true,cppqedutils::trajectory::observerNoOp);*/

  // ****** Parameters of the Problem
  
  /*ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  ParsDrivenDissipative pplqb(p); */

  // Parameter finalization
  // QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  /*
  StateVector psi(state0()()+state1()());
  psi.renorm();
  */
  
  // evolve(init(pplqb),make(pplqb,qmp),pe);

}



/*
  static_assert( hamiltonian<decltype(diagonalH(z)),1> ) ;
  static_assert( system_frequency_store< decltype(std::array{structure::SystemFrequencyDescriptor{"z",z,1.},structure::SystemFrequencyDescriptor{"η",eta,1.}}) > );
  static_assert( expectationvalues::time_independent_functional<decltype(expectationValues),1> ) ;
*/
