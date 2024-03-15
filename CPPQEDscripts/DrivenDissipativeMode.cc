// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// #include "Evolution.h"
#include "Mode.h"

#include "Master.h"
#include "QuantumJumpMonteCarlo.h"

#include <iostream>


using namespace mode;
using namespace quantumtrajectory;

static const auto algorithm=qjmc::Algorithm::stepwise;

int main(int argc, char* argv[])
{
  auto op{optionParser()};

  trajectory::Pars<qjmc::Pars<algorithm,pcg64>> pt(op);

  Pars pm(op);

  parse(op,argc, argv);

  auto mode{make(pm)};

  quantumdata::StateVector<1> psi{{pm.cutoff}}; psi(0)=0; psi(9)=1;

  run<algorithm,ODE_EngineBoost>(std::move(mode),std::move(psi),pt,trajectory::observerNoOp);

  // run(
  //   // qjmc::make<algorithm,ODE_EngineBoost>(std::move(mode),std::move(psi),pt),
  //   qjmc::makeEnsemble<algorithm,ODE_EngineBoost>(std::move(mode),psi,pt),
  //   pt,trajectory::observerNoOp);

 // std::cerr<<json(nOp(10))<<std::endl;
}

