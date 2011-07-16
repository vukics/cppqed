#include "BichromaticMode.h"

#include "Pars.h"

mode::ParsBichromatic::ParsBichromatic(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod), ParsPumpedLossy(p,mod),
    deltaOther(p.addTitle("BichromaticMode",mod).addMod("deltaC_Other",mod,"Other mode detuning",0.)),
    etaOther  (p.addMod("etaOther",mod,"Other cavity pump",dcomp(0)))
{}

