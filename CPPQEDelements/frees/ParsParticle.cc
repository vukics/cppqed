// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsParticle.h"

#include "ParticleInitialCondition.h"

#include "Pars.tcc"


namespace particle {

Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : omrec(p.addTitle("Particle",mod).addMod("omrec",mod,"Particle recoil frequency",1.)),
    fin(p.addMod<size_t>("fin",mod,"Particle space resolution: 1<<fin",6)),
    init(p.addMod("pinit",mod,"Wavepacket initial condition for particle",InitialCondition(-.5,0,.1))),
    hoInitn(p.addMod<int>("HOinitn",mod,"Particle initial condition harmonic oscillator nth eigenstate with char. length pinit.sig",-1))
{
}


ParsPumped::ParsPumped(parameters::ParameterTable& p, const std::string& mod)
  : Pars(p,mod),
    vClass(p.addTitle("PumpedParticle",mod).addMod("vClass",mod,"Particle effective pumping strength",-10.)),
    kPart(p.addMod<size_t>("kPart",mod,"Particle pump wavenumber",1)),
    modePart(p.addMod("modePart",mod,"Particle pump mode function",MFT_SIN))
{
}

} // particle
