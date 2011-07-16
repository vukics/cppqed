#include "ParsMode.h"

#include "Pars.h"


namespace mode {


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : cutoff(p.addTitle("Mode",mod).addMod<size_t>("cutoff",mod,"Mode cutoff",10)),
    minitFock(p.addMod<size_t>("minitFock",mod,"Mode initial Fock state",0)),
    minit(p.addMod<dcomp>("minit",mod,"Mode initial field",0)),
    displayLevel(p.addMod<size_t>("modeDisplayLevel",mod,"When nonzero, quadrature variances are also calculated",0)),
    delta(p.addMod("deltaC",mod,"Mode detuning",-10.))
{}


ParsPumped::ParsPumped(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod),
    eta(p.addTitle("PumpedMode",mod).addMod("eta",mod,"Cavity pump",dcomp(0)))
{}


ParsLossy::ParsLossy(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod),
    kappa(p.addTitle("LossyMode",mod).addMod("kappa",mod,"Mode decay rate",-delta)),
    nTh  (p.addMod("nTh"  ,mod,"Mode thermal-photon number",0.))
{}


ParsPumpedLossy::ParsPumpedLossy(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod), ParsPumped(p,mod), ParsLossy(p,mod)
{}


} // mode
