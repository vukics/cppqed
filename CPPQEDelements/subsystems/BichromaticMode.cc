// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BichromaticMode.h"

#include "Pars.h"

mode::ParsBichromatic::ParsBichromatic(parameters::Table& p, const std::string& mod)
  : Pars(p,mod), ParsDrivenDissipative(p,mod),
    deltaOther(p.addTitle("BichromaticMode",mod).add("deltaC_Other",mod,"Other mode detuning",0.)),
    etaOther  (p.add("etaOther",mod,"Other cavity pump",dcomp(0)))
{}

