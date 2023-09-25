// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// #include "Evolution.h"
#include "Mode.h"

#include "Master.h"
#include "QuantumJumpMonteCarlo.h"

#include <iostream>


using cppqedutils::json;
using namespace mode;

int main(int argc, char* argv[])
{
  auto op{optionParser()};

  cppqedutils::trajectory::Pars<quantumtrajectory::qjmc::Pars<pcg64>> pt(op);

  Pars pm(op);

  parse(op,argc, argv);

  auto mode{make(pm)};

  quantumdata::StateVector<1> psi{{pm.cutoff}}; psi(0)=0; psi(1)=1;

  cppqedutils::run(
    // quantumtrajectory::qjmc::make<cppqedutils::ODE_EngineBoost>(std::move(mode),/*quantumdata::StateVector<1>{{pm.cutoff}}*/std::move(psi),pt),
    quantumtrajectory::qjmc::makeEnsemble<cppqedutils::ODE_EngineBoost>(std::move(mode),/*quantumdata::StateVector<1>{{pm.cutoff}}*/std::move(psi),pt),
    pt,cppqedutils::trajectory::observerNoOp);
}

