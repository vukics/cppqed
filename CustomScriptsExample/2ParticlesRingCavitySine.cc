// Copyright Raimar Sandner 2012–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"
#include "Mode.h"
#include "Particle.h"
#include "ParticleCavity.h"
#include "ParticleCavity_InterferenceWorkaround.h"
#include "MomentumCorrelation.h"


using namespace quantumdata;
using namespace particlecavity_interferenceworkaround;

int main(int argc, char **argv)
{
    ParameterTable p;

    evolution::Pars<> pe(p);              // Driver parameters
    mode::ParsDissipative plm(p);             // Dissipative mode (cosine)
    particle::ParsDriven ppp1(p,"1");   // Driven particle 1 (sine)
    particle::ParsDriven ppp2(p,"2");   // Driven particle 2
    particlecavity::ParsAlong ppci(p);  // Particle-Cavity interaction

    bool &bosons = p.add("bosons", "Particle are bosons", false);
    bool &fermions = p.add("fermions", "Particles are fermions", false);
    
    double &alphasquare = p.add("alphasquare", "Coherent state of driven cavity mode (has to be real)", 1.);
    pe.evol = evolution::SINGLE;
    
    QM_Picture& qmp=updateWithPicture(p,argc,argv);
    
    // Parameters, take U0 as configuration parameter and update V0 and Delta accordingly
    
    ppci.modeCav = MFT_COS;
    ppp1.vClass = ppp2.vClass = ppci.uNot * alphasquare;
    plm.delta -= ppci.uNot;

    particle::PtrDriven myParticle(particle::makeDriven(ppp1,qmp));
    mode::Ptr myMode(mode::make(plm,qmp));
    
    auto particlecavity{std::make_shared<ParticleAlongCavity>(myMode,myParticle,ppci)};
    auto mci{std::make_shared<MomentumCorrelation>(myParticle,myParticle)};
    
    StateVector<1> stateParticle1(particle::init(ppp1));
    StateVector<1> stateParticle2(particle::init(ppp2));
    StateVector<2> stateParticles = stateParticle1*stateParticle2;
    
    if (bosons) {
        stateParticles = stateParticle1*stateParticle2 + stateParticle2*stateParticle1;
    }
    else if (fermions) {
        stateParticles = stateParticle1*stateParticle2 - stateParticle2*stateParticle1;
    }
    stateParticles.renorm();
    
    StateVector<3> psi(mode::init(plm)*stateParticles);
    
    evolve(psi,composite::make(_<0,1>(particlecavity),_<0,2>(particlecavity),
                               _<1,2>(mci)),pe);


}
