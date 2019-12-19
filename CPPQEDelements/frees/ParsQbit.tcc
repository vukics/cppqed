// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsQbit.h"

#include "Pars.tcc"


namespace qbit {


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : qbitInit(p.addTitle("Qbit",mod).addMod("qbitInit",mod,"Qbit initial condition excited",dcomp(0.))),
    delta(p.addMod("deltaA",mod,"Qbit detuning",-10.))
{}


ParsPumped::ParsPumped(parameters::ParameterTable& p, const std::string& mod)
    : Pars(p,mod), eta(p.addTitle("PumpedQbit",mod).addMod("etat",mod,"Qbit pump",dcomp(0.))) {}


template <typename BASE>
ParsLossy<BASE>::ParsLossy(parameters::ParameterTable& p, const std::string& mod)
    : BASE(p,mod), gamma(p.addTitle("LossyQbit",mod).addMod("gamma",mod,"Qbit decay rate",fabs(BASE::delta))) {}


template <typename BASE>
ParsLossyPhaseNoise<BASE>::ParsLossyPhaseNoise(parameters::ParameterTable& p, const std::string& mod)
    : BASE(p,mod), gamma_parallel(p.addMod("gamma_parallel",mod,"Qbit phase flip rate",0.)) {}


} // qbit
