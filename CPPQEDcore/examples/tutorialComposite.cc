#include "EvolutionComposite.h"

#include "ParticleTwoModes.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsEvolution pe(p);
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

  evolve(psi,
         composite::make(         
			 Act<1,0>  (ParticleAlongCavity(plus ,part,ppcP)),
			 Act<2,0>  (ParticleAlongCavity(minus,part,ppcM)),
			 Act<1,2,0>(ParticleTwoModes(plus,minus,part,ppcP,ppcM))
			 ),
         pe);

}