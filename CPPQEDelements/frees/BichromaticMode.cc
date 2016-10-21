// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BichromaticMode.h"

#include "Pars.tcc"

mode::ParsBichromatic::ParsBichromatic(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod), ParsPumpedLossy(p,mod),
    deltaOther(p.addTitle("BichromaticMode",mod).addMod("deltaC_Other",mod,"Other mode detuning",0.)),
    etaOther  (p.addMod("etaOther",mod,"Other cavity pump",dcomp(0)))
{}

