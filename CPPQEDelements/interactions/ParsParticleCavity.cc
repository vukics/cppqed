// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsParticleCavity.h"

#include "Pars.h"


particlecavity::ParsGenericPump::ParsGenericPump(parameters::Table& p, const std::string& mod)
  : etaeff(p.addMod("etaeff",mod,"Particle-Cavity interaction effective cavity pump parameter",0.)) {}


particlecavity::ParsOrthogonal ::ParsOrthogonal(parameters::Table& p, const std::string& mod)
  : uNot(p.addTitle("ParticleCavityOrthogonal",mod).addMod("Unot",mod,"Particle-Cavity interaction Unot parameter",-1.)) {}



particlecavity::ParsAlong      ::ParsAlong     (parameters::Table& p, const std::string& mod)
  : ParsOrthogonal(p,mod),
    kCav(p.addTitle("ParticleCavityAlong",mod).addMod<size_t>("kCav",mod,"Cavity mode wavenumber",1)),
    modeCav(p.addMod("modeCav",mod,"Cavity mode function",MFT_SIN)) {}
