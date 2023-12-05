// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// #include "Evolution.h"
#include "Spin.h"

using namespace std ; using cppqedutils::json ;


int main(int argc, char* argv[])
{
//   // ****** Parameters of the Problem
//
//   ParameterTable p;
//
//   evolution::Pars<> pe(p); // Driver Parameters
//
//   spin::Pars ps(p);
//
//   // Parameter finalization
//   update(p,argc,argv,"--");
//
//   // ****** ****** ****** ****** ****** ******
//
//   const auto spin{std::make_shared<DissipativeSpinSch>(ps)};
//
//   structure::freesystem::StateVector psi{spin->getDimensions()};
//
//   psi(psi.getArray().ubound(0))=1;
//
//   evolve(psi,spin,pe);

//  auto sp=spin::splus(2)//, sz=spin::sz(2) ;

  std::cerr<<json(spin::splus(2))<<std::endl<<json(spin::sminus(2))<<std::endl<<json(spin::sz(2))<<std::endl;

}

