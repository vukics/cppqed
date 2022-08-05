// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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


struct ParsPumped : Pars
{
  dcomp& eta;

  ParsPumped(parameters::Table& p, const std::string& mod="")
    : Pars(p,mod), eta(p.addTitle("PumpedQbit",mod).add("etat",mod,"Qbit pump",dcomp(0.))) {}

};


template <typename BASE=Pars>
struct ParsLossy : BASE
{
  double &gamma;

  ParsLossy(parameters::Table& p, const std::string& mod="")
    : BASE(p,mod), gamma(p.addTitle("LossyQbit",mod).add("gamma",mod,"Qbit decay rate",fabs(BASE::delta))) {}
};


typedef ParsLossy<ParsPumped> ParsPumpedLossy;


template <typename BASE=ParsLossy<>>
struct ParsLossyPhaseNoise : BASE
{
  double &gamma_parallel;
  
  ParsLossyPhaseNoise(parameters::Table& p, const std::string& mod="")
    : BASE(p,mod), gamma_parallel(p.add("gamma_parallel",mod,"Qbit phase flip rate",0.)) {}
  
};


typedef ParsLossyPhaseNoise<ParsPumpedLossy> ParsPumpedLossyPhaseNoise;


template <typename BASE>
struct ParsLossyIncoherentPump : BASE
{
  double &gamma_pump;
  
  ParsLossyIncoherentPump(parameters::Table& p, const std::string& mod="")
    : BASE(p,mod), gamma_pump(p.add("gamma_pump",mod,"Qbit incoherent pump rate",0.)) {}
  
};


} // qbit

#endif // CPPQEDELEMENTS_FREES_PARSQBIT_H_INCLUDED
