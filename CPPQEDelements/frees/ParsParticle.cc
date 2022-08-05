// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsParticle.h"

#include "ParticleInitialCondition.h"

#include "Pars.h"


namespace particle {

Pars::Pars(parameters::Table& p, const std::string& mod)
  : omrec(p.addTitle("Particle",mod).add("omrec",mod,"Particle recoil frequency",1.)),
    fin(p.add<size_t>("fin",mod,"Particle space resolution: 1<<fin",6)),
    init(p.add("pinit",mod,"Wavepacket initial condition for particle",InitialCondition(-.5,0,.1))),
    hoInitn(p.add<int>("HOinitn",mod,"Particle initial condition harmonic oscillator nth eigenstate with char. length pinit.sig",-1))
{
}


ParsDriven::ParsDriven(parameters::Table& p, const std::string& mod)
  : Pars(p,mod),
    vClass(p.addTitle("DrivenParticle",mod).add("vClass",mod,"Particle effective pumping strength",-10.)),
    kPart(p.add<size_t>("kPart",mod,"Particle pump wavenumber",1)),
    modePart(p.add("modePart",mod,"Particle pump mode function",MFT_SIN))
{
}

} // particle
