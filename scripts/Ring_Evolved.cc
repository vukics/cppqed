#include "SimulatedHigh.h"

#include "ParticleTwoModes.h"

#include "ParticleInitialCondition.h"


using namespace std       ;
using namespace cpputils  ;
using namespace trajectory;
using namespace mathutils ;

typedef TTD_CARRAY(1) Array;


void derivs(double, const Array& b, Array& dbdt,
	    const particle::Pars& pp,
	    const mode::ParsPumpedLossy& pm1, const mode::ParsPumpedLossy& pm2,
	    const particlecavity::ParsAlong pac1, const particlecavity::ParsAlong pac2)
{
  static const double x=PI*pp.init.getX0();

  static const dcomp z1=DCOMP_I*(pm1.delta-pac1.uNot*sqrAbs(modeFunction(pac1.modeCav,x)))-pm1.kappa;
  static const dcomp z2=DCOMP_I*(pm2.delta-pac2.uNot*sqrAbs(modeFunction(pac2.modeCav,x)))-pm2.kappa;

  static const dcomp g=sign(pac1.uNot)*sqrt(pac1.uNot*pac2.uNot)*
    conj(modeFunction(pac1.modeCav,x))*modeFunction(pac2.modeCav,x);

  dbdt(0)=z1*b(0)+pm1.eta-DCOMP_I*     g *b(1);
  dbdt(1)=z2*b(1)+pm2.eta-DCOMP_I*conj(g)*b(0);

}


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsTrajectory pt(p);

  particle::Pars pp(p);
  mode::ParsPumpedLossy pmP(p,"P");
  mode::ParsPumpedLossy pmM(p,"M");
  particlecavity::ParsAlong ppcP(p,"P");
  particlecavity::ParsAlong ppcM(p,"M");

  ppcP.modeCav=MFT_PLUS; ppcM.modeCav=MFT_MINUS; 

  update(p,argc,argv,"--");

  Array alpha(2);
  
  alpha(0)=pmP.minit;
  alpha(1)=pmM.minit;

  Simulated<Array> S(alpha,bind(derivs,_1,_2,_3,pp,pmP,pmM,ppcP,ppcM),1e-6,Array(),pt);

  evolve(S,pt);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
