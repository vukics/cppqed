#include "ParsQbit.h"

#include "impl/Pars.tcc"


namespace qbit {


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : qbitInit(p.addTitle("Qbit",mod).addMod("qbitInit",mod,"Qbit initial condition excited",dcomp(0.))),
    delta(p.addMod("deltaA",mod,"Qbit detuning",-10.))
{}


ParsPumped::ParsPumped(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod),
    eta(p.addTitle("PumpedQbit",mod).addMod("etat",mod,"Qbit pump",dcomp(0)))
{}


ParsLossy::ParsLossy(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod),
    gamma(p.addTitle("LossyQbit",mod).addMod("gamma",mod,"Qbit decay rate",-delta))
{}


ParsPumpedLossy::ParsPumpedLossy(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod), ParsPumped(p,mod), ParsLossy(p,mod)
{}


} // qbit
