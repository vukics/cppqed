// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsParticleCavity.h"

#include "Pars.tcc"


namespace particlecavity {

ParsOrthogonal::ParsOrthogonal(parameters::ParameterTable& p, const std::string& mod)
  : uNot(p.addTitle("ParticleCavityOrthogonal",mod).addMod("Unot",mod,"Particle-Cavity interaction Unot parameter",-1.)) {}


ParsAlong     ::ParsAlong     (parameters::ParameterTable& p, const std::string& mod)
  : ParsOrthogonal(p,mod),
    kCav(p.addTitle("ParticleCavityAlong",mod).addMod<size_t>("kCav",mod,"Cavity mode wavenumber",1)),
    modeCav(p.addMod("modeCav",mod,"Cavity mode function",MFT_SIN)) {}


} // particlecavity
