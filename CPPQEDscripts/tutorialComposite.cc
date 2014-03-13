#include "EvolutionComposite.h"

#include "ParticleTwoModes.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars pe(p);
  particle::Pars pp(p);
  mode::ParsLossy       pmP(p,"P");
  mode::ParsPumpedLossy pmM(p,"M");
  particlecavity::ParsAlong ppcP(p,"P");
  particlecavity::ParsAlong ppcM(p,"M");

  ppcP.modeCav=MFT_PLUS; ppcM.modeCav=MFT_MINUS; 

  update(p,argc,argv,"--");

  particle::Ptr part (make(pp ,QMP_IP));
  mode    ::Ptr plus (make(pmP,QMP_IP));
  mode    ::Ptr minus(make(pmM,QMP_IP));

  quantumdata::StateVector<3> psi(wavePacket(pp)*init(pmP)*init(pmM));

  evolve<0>(psi,
            composite::make(
                            _<1,0>  (ParticleAlongCavity(plus ,part,ppcP)),
                            _<2,0>  (ParticleAlongCavity(minus,part,ppcM)),
                            _<1,2,0>(ParticleTwoModes(plus,minus,part,ppcP,ppcM))
                            ),
            pe);

  {
  ParticleAlongCavity pac1(plus ,part,ppcP), pac2(minus,part,ppcM);
  ParticleTwoModes ptm(plus,minus,part,ppcP,ppcM);

  const composite::result_of::Make<_<1,0>,_<2,0>,_<1,2,0> >::type
    system=composite::make(
                           _<1,0>  (pac1),
                           _<2,0>  (pac2),
                           _<1,2,0>(ptm));

  system->displayParameters(std::cout);
}
  {
  ParticleAlongCavity pac1(plus ,part,ppcP), pac2(minus,part,ppcM);
  ParticleTwoModes ptm(plus,minus,part,ppcP,ppcM);

  const auto
    system=composite::make(
                           _<1,0>  (pac1),
                           _<2,0>  (pac2),
                           _<1,2,0>(ptm));

  system->displayParameters(std::cout);
}

}


#include "IdenticalParticles.h"

void f(int argc, char* argv[])
{
ParameterTable p;

evolution::Pars pe(p);
particle::ParsPumped pp(p);
mode::ParsLossy pm(p);
particlecavity::ParsOrthogonal ppc(p);

update(p,argc,argv,"--");

LossyMode<>    mode(pm); // Free0
PumpedParticle part(pp); // Free1,2 – only one instant

ParticleOrthogonalToCavity act(mode,part,ppc); // only one instant

quantumdata::StateVector<3> psi(init(pm)*hoState(pp)*hoState(pp));

evolve(psi,
       composite::make(
                       _<0,1>(act),
                       _<0,2>(act),
                       _<1,2>(IdenticalParticles<2>(part /*,...*/))
                       ),
       pe);
}
