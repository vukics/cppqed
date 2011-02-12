#include "EvolutionHigh.h"

#include "GeneralDicke.h"

#include "BinarySystem.h"


using namespace std;

typedef quantumdata::StateVector<1> StateVector1;
typedef quantumdata::StateVector<2> StateVector2;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters

  mode::ParsPumpedLossy pplm(p);

  spin::Pars ps(p);

  // generaldicke::Pars pgd(p);

  double& u=p.add("u","General Dicke interaction u parameter",1.);
  double& y=p.add("y","General Dicke interaction y parameter",1.);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  u/=ps.twos; y/=sqrt(ps.twos);

  mode::SmartPtr mode(mode::maker(pplm,QMP_IP));

  Spin spin(ps);

  GeneralDicke<> gd(mode,spin,u,y);

  BinarySystem sys(gd);

  StateVector1 psiMode(mode::init(pplm)), psiSpin(spin.getDimension());

  psiSpin()(0)=1; psiSpin.renorm();

  StateVector2 psi(psiMode*psiSpin);

  evolve(psi,sys,pe,tmptools::Vector<0>());

  } catch (ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
