// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// #include "Evolution.h"

#include "Qbit.h"

#include <iostream>

using namespace qbit;
using namespace std;

void f(hamiltonian<1> auto) {}

int main(int argc, char* argv[])
{
  dcomp z{.1,1.}, eta{-.1,.2};

  static_assert( hamiltonian<decltype(diagonalH(z)),1> ) ;
  static_assert( system_frequency_store< decltype(std::array{structure::SystemFrequencyDescriptor{"z",z,1.},structure::SystemFrequencyDescriptor{"η",eta,1.}}) > );
  static_assert( expectationvalues::time_independent_functional<decltype(expectationValues),1> ) ;

  QuantumSystemDynamics qbit {
    std::array{structure::SystemFrequencyDescriptor{"z",z,1.},structure::SystemFrequencyDescriptor{"η",eta,1.}},
    std::array{loss(1.),gain(.1)},
    diagonalH(z),
    expectationValues }; //= make({1,1},{1,1},{1,1},1,1,1);
  
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

