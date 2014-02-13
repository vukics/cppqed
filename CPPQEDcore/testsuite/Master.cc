#include "Evolution.h"

#include "Mode.h"
#include "Qbit.h"

#include "JaynesCummings.h"

#include "BinarySystem.h"

using namespace std;
using namespace mode;
using namespace quantumtrajectory;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  evolution::Pars pe(p); // Driver Parameters
  ParsPumpedLossy pplm(p); 

  qbit::ParsPumpedLossy pplqb(p); 
  jaynescummings::Pars  pjc  (p); 

  pplm.kappa=0;

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if (qmp==QMP_IP) qmp=QMP_UIP;

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  mode::Ptr mode(mode::make(pplm ,qmp));

  JaynesCummings<> jc(qbit,mode,pjc);

  BinarySystem system(jc);

  {
    StateVector     psi(mode::coherent(dcomp(1,-.5),pplm.cutoff)), dpsidt(psi.getDimensions());
    DensityOperator rho(psi), drhodt(rho.getDimensions());

    psi.renorm(); rho.renorm();

    MCWF_Trajectory<1> trajS(psi,*mode,pe);
    Master<1>          trajM(rho,*mode,pe,false);

    trajS.derivs(pe.T,psi(),dpsidt()); psi()+=pe.Dt*dpsidt();
    trajM.derivs(pe.T,rho(),drhodt()); rho()+=pe.Dt*drhodt();
    
    trajS.derivs(pe.T,psi(),dpsidt()); psi()+=pe.Dt*dpsidt();
    trajM.derivs(pe.T,rho(),drhodt()); rho()+=pe.Dt*drhodt();
    
    rho()-=psi.dyad();
    // cout<<drhodt()<<dpsidt.dyad();
    cout<<max(abs(rho()))<<endl;
  }

  {
    quantumdata::StateVector    <2> psi(qbit::init(pplqb)*mode::init(pplm)), dpsidt(psi.getDimensions());
    quantumdata::DensityOperator<2> rho(psi), drhodt(rho.getDimensions());

    psi.renorm(); rho.renorm();

    MCWF_Trajectory<2> trajS(psi,system,pe);
    Master<2>          trajM(rho,system,pe,false);

    trajS.derivs(pe.T,psi(),dpsidt()); psi()+=pe.Dt*dpsidt();
    trajM.derivs(pe.T,rho(),drhodt()); rho()+=pe.Dt*drhodt();

    rho()-=psi.dyad();
    // cout<<drhodt()<<dpsidt.dyad();
    cout<<max(abs(rho()))<<endl;

    pe.T=0;

    evolve(psi,system,pe);

  }

  // StateVector psi(mode::fock(2,pplm.cutoff));

  /*
  cout<<psi();
  structure::Liouvillean<1>::actWithJ(&plm,psi(),0);
  psi()/=dcomp(1,-.5); psi.renorm();
  cout<<psi();
  */


  } catch (NamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
