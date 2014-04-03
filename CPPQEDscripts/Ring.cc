// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"

#include "ParticleTwoModes.h"


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars pe(p);
  particle::Pars pp(p);
  mode::ParsPumpedLossy pmP(p,"P");
  mode::ParsPumpedLossy pmM(p,"M");
  particlecavity::ParsAlong ppcP(p,"P");
  particlecavity::ParsAlong ppcM(p,"M");

  ppcP.modeCav=MFT_PLUS; ppcM.modeCav=MFT_MINUS; 

  QM_Picture& qmp=updateWithPicture(p,argc,argv);

  pmP.delta-=ppcP.uNot/(isComplex(ppcP.modeCav) ? 1. : 2.);
  pmM.delta-=ppcM.uNot/(isComplex(ppcM.modeCav) ? 1. : 2.);

  particle::Ptr part (make(pp ,qmp));
  mode    ::Ptr plus (make(pmP,qmp));
  mode    ::Ptr minus(make(pmM,qmp));

  quantumdata::StateVector<3> psi(init(pp)*
                                  init(pmP)*
                                  init(pmM));

  ParticleAlongCavity pacP(plus ,part,ppcP);
  ParticleAlongCavity pacM(minus,part,ppcM);
  ParticleTwoModes ptm(plus,minus,part,ppcP,ppcM);

  evolve<0>(psi,
            composite::make(
                            _<1,0>  (pacP),
                            _<2,0>  (pacM),
                            _<1,2,0>(ptm)
                            ),
            pe);

}
