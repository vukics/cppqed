// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsMode.h"

#include "Pars.h"


namespace mode {


Pars::Pars(parameters::Table& p, const std::string& mod)
  : cutoff(p.addTitle("Mode",mod).add<size_t>("cutoff",mod,"Mode cutoff",10)),
    minitFock(p.add<size_t>("minitFock",mod,"Mode initial Fock state",0)),
    minit(p.add<dcomp>("minit",mod,"Mode initial field",0)),
    streamLevel(p.add<size_t>("modeStreamLevel",mod,"When nonzero, quadrature variances are also calculated",0)),
    delta(p.add("deltaC",mod,"Mode detuning",-10.)),
    omegaKerr(p.addTitle("Kerr nonlinearity",mod).add("omegaKerr",mod,"Kerr constant",0.)),
    omegaKerrAlter(p.add("omegaKerrAlter",mod,"Kerr constant – alternative implementation",0.))
{}


ParsPumped::ParsPumped(parameters::Table& p, const std::string& mod)
  : Pars(p,mod),
    eta(p.addTitle("PumpedMode",mod).add("eta",mod,"Cavity pump",dcomp(0)))
{}


ParsLossy::ParsLossy(parameters::Table& p, const std::string& mod)
  : Pars(p,mod),
    kappa(p.addTitle("LossyMode",mod).add("kappa",mod,"Mode decay rate",-delta)),
    nTh  (p.add("nTh"  ,mod,"Mode thermal-photon number",0.))
{}


ParsPumpedLossy::ParsPumpedLossy(parameters::Table& p, const std::string& mod)
  : Pars(p,mod), ParsPumped(p,mod), ParsLossy(p,mod)
{}


} // mode
