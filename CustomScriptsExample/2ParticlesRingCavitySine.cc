#include "EvolutionComposite.h"
#include "Mode.h"
#include "Particle.h"
#include "ParticleCavity.h"
#include "ParticleCavity_InterferenceWorkaround.h"
#include "MomentumCorrelation.h"
#include "Version.h"

#include "component_versions.h"

using namespace quantumdata;
using namespace particlecavity_interferenceworkaround;

int main(int argc, char **argv)
{
    updateVersionstring(cppqed_component_versions());
    ParameterTable p;

    evolution::Pars pe(p);              // Driver parameters
    mode::ParsLossy plm(p);             // Lossy mode (cosine)
    particle::ParsPumped ppp1(p,"1");   // Pumped particle 1 (sine)
    particle::ParsPumped ppp2(p,"2");   // Pumped particle 2
    particlecavity::ParsAlong ppci(p);  // Particle-Cavity interaction
    ParsInterference pi(p);             // Interference term
    
    bool &bosons = p.add("bosons", "Particle are bosons", false);
    bool &fermions = p.add("fermions", "Particles are fermions", false);
    
    double &alphasquare = p.add("alphasquare", "Coherent state of pumped cavity mode (has to be real)", 1.);
    QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);
    pe.evol = evolution::SINGLE;
    
    update(p,argc,argv,"--");
    if (pe.evol == evolution::MASTER && qmp != QMP_UIP) qmp = QMP_UIP;
    
    // Parameters, take U0 as configuration parameter and update V0 and Delta accordingly
    
    ppp1.vClass = ppp2.vClass = ppci.uNot * alphasquare;
    plm.delta -= ppci.uNot;
    pi.uInterference = ppci.uNot*sqrt(alphasquare)/2.;
    
    particle::Ptr myParticle(particle::make(ppp1,qmp));
    mode::Ptr myMode(mode::make(plm,qmp));
    
    ParticleAlongCavity particlecavity(myMode,myParticle,ppci,0);
    Interference interference(myMode,myParticle,pi);
    MomentumCorrelation mci(myParticle,myParticle);
    
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
    
    evolve(psi,composite::make(Act<0,1>(particlecavity),Act<0,2>(particlecavity),
                               Act<0,1>(interference),Act<0,2>(interference),Act<1,2>(mci)),pe);
}
