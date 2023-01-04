// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsParticleCavity.h"

#include "Pars.h"


particlecavity::ParsGenericPump::ParsGenericPump(parameters::Table& p, const std::string& mod)
  : etaeff(p.add("etaeff",mod,"Particle-Cavity interaction effective cavity pump parameter",0.)) {}


particlecavity::ParsOrthogonal ::ParsOrthogonal(parameters::Table& p, const std::string& mod)
  : uNot(p.addTitle("ParticleCavityOrthogonal",mod).add("Unot",mod,"Particle-Cavity interaction Unot parameter",-1.)) {}



particlecavity::ParsAlong      ::ParsAlong     (parameters::Table& p, const std::string& mod)
  : ParsOrthogonal(p,mod),
    kCav(p.addTitle("ParticleCavityAlong",mod).add<size_t>("kCav",mod,"Cavity mode wavenumber",1)),
    modeCav(p.add("modeCav",mod,"Cavity mode function",MFT_SIN)) {}
