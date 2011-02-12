#include "EvolutionHigh.h"
#include "Composite.h"
#include "ParticleCavity.h"
#include "ParticleTwoModes.h"

#include "ModeFunctionType.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsEvolution pe(p);
  particle::Pars pp(p);
  mode::ParsPumpedLossy pmP(p,"P");
  mode::ParsPumpedLossy pmM(p,"M");
  particlecavity::ParsAlong ppcP(p,"P");
  particlecavity::ParsAlong ppcM(p,"M");

  ppcP.modeCav=MFT_PLUS; ppcM.modeCav=MFT_MINUS; 

  update(p,argc,argv,"--");

  pmP.delta-=ppcP.uNot/(isComplex(ppcP.modeCav) ? 1. : 2.);
  pmM.delta-=ppcM.uNot/(isComplex(ppcM.modeCav) ? 1. : 2.);

  particle::SmartPtr part (maker(pp ,QMP_IP));
  mode    ::SmartPtr plus (maker(pmP,QMP_IP));
  mode    ::SmartPtr minus(maker(pmM,QMP_IP));

  quantumdata::StateVector<3> psi(init(pp)*
				  init(pmP)*
				  init(pmM));

  ParticleAlongCavity pacP(plus ,part,ppcP);
  ParticleAlongCavity pacM(minus,part,ppcM);
  ParticleTwoModes ptm(plus,minus,part,ppcP,ppcM);

  evolve(psi,
	 makeComposite(	         
		       Act<1,0>  (pacP),
		       Act<2,0>  (pacM),
		       Act<1,2,0>(ptm)
				 ),
	 pe);

}


/*
  {

    quantumdata::StateVector<3> psi(init(pmP)*init(pmM)*init(pp));

    psi()=0;
  
    quantumdata::StateVector<3> dpsidt(psi.getDimensions());
  
    psi()(0,1,7)=1;

    static_cast<structure::Hamiltonian<3>*>(&ptm)->addContribution(0,psi(),dpsidt(),0);

    std::cout<<dpsidt();
  
  }
*/
