// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// #include "Evolution.h"

#include "Qbit.h"

#include "Master.h"
#include "QuantumJumpMonteCarlo.h"

#include <iostream>

using namespace qbit;
using namespace std;



int main(int argc, char* argv[])
{
  auto op{optionParser()};

  cppqedutils::trajectory::Pars<quantumtrajectory::qjmc::Pars<pcg64>> pt(op);

  Pars pq(op);

  parse(op,argc, argv);

  auto qbit{make(pq)};

  cppqedutils::ODE_EngineBoost<quantumdata::DensityOperator<1>::StorageType> oe(quantumtrajectory::initialTimeStep(getFreqs(qbit)),1e-12,1e-30);

  quantumtrajectory::Master m{qbit,quantumdata::DensityOperator<1>{{2}},oe};

  // quantumtrajectory::QuantumJumpMonteCarlo q{qbit,quantumdata::StateVector<1>{{2}},oe,randomutils::EngineWithParameters<pcg64>{1001,1},0.01};

  cppqedutils::run(m,pt,cppqedutils::trajectory::observerNoOp);

  // {
  //   quantumdata::StateVector<2> psi{{2,3}}; quantumdata::DensityOperator<2> rho{{2,3}};
  //   std::cerr<<_(psi,{0,1},{1,0})<<" "<<_(rho,{0,1},{1,0})<<std::endl;
  // }


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
