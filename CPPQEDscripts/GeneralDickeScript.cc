#include "EvolutionBinary.h"
#include "Mode.h"

#include "GeneralDicke.h"


using namespace std;

typedef quantumdata::StateVector<1> StateVector1;
typedef quantumdata::StateVector<2> StateVector2;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters

  mode::ParsPumpedLossy pplm(p);

  spin::Pars ps(p);

  // generaldicke::Pars pgd(p);

  double& u=p.addTitle("Script").add("u","General Dicke interaction u parameter",1.);
  double& y=p.add("y","General Dicke interaction y parameter",1.);

  QM_Picture& qmp=p.add("picture","QM_Picture for mode (IP=UIP or Sch) and field (IP, UIP, or Sch)",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  u/=ps.twoS; y/=sqrt(ps.twoS);

  Spin spin(ps);

  StateVector1 psiMode(mode::init(pplm)), psiSpin(spin.getDimension());

  psiSpin.getArray()(0)=1; psiSpin.renorm();

  StateVector2 psi(psiMode*psiSpin);

  evolve<0>
    (psi,
     binary::make(GeneralDicke<>(mode::make<mode::AveragedMonitorCutoff<mode::AveragedQuadratures> >(pplm,qmp),spin,u,y)),
     pe);




}
