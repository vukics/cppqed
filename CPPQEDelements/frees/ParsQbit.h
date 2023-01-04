// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_PARSQBIT_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSQBIT_H_INCLUDED

#include "Pars.h"

#include "ComplexExtensions.h"


namespace qbit {

struct Pars
{
  dcomp &qbitInit;
  double &delta;

  Pars(parameters::Table& p, const std::string& mod="")
    : qbitInit(p.addTitle("Qbit",mod).add("qbitInit",mod,"Qbit initial condition excited",dcomp(0.))),
      delta(p.add("deltaA",mod,"Qbit detuning",-10.)) {}
      
};


struct ParsDriven : Pars
{
  dcomp& eta;

  ParsDriven(parameters::Table& p, const std::string& mod="")
    : Pars(p,mod), eta(p.addTitle("DrivenQbit",mod).add("etat",mod,"Qbit pump",dcomp(0.))) {}

};


template <typename BASE=Pars>
struct ParsDissipative : BASE
{
  double &gamma;

  ParsDissipative(parameters::Table& p, const std::string& mod="")
    : BASE(p,mod), gamma(p.addTitle("DissipativeQbit",mod).add("gamma",mod,"Qbit decay rate",fabs(BASE::delta))) {}
};


typedef ParsDissipative<ParsDriven> ParsDrivenDissipative;


template <typename BASE=ParsDissipative<>>
struct ParsDissipativePhaseNoise : BASE
{
  double &gamma_parallel;
  
  ParsDissipativePhaseNoise(parameters::Table& p, const std::string& mod="")
    : BASE(p,mod), gamma_parallel(p.add("gamma_parallel",mod,"Qbit phase flip rate",0.)) {}
  
};


typedef ParsDissipativePhaseNoise<ParsDrivenDissipative> ParsDrivenDissipativePhaseNoise;


template <typename BASE>
struct ParsDissipativeIncoherentPump : BASE
{
  double &gamma_pump;
  
  ParsDissipativeIncoherentPump(parameters::Table& p, const std::string& mod="")
    : BASE(p,mod), gamma_pump(p.add("gamma_pump",mod,"Qbit incoherent pump rate",0.)) {}
  
};


} // qbit

#endif // CPPQEDELEMENTS_FREES_PARSQBIT_H_INCLUDED
