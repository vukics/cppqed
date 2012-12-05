#include "EvolutionBinary.h"
#include "Mode.h"
#include "Spin.h"

#include "JaynesCummings.h"

#include "MatrixOfHamiltonian.h"


using namespace std;

typedef quantumdata::StateVector<2> StateVector2;
typedef quantumdata::StateVector<1> StateVector1;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
//  qbit::ParsPumpedLossy pplqb(p); 

  spin::Pars sp (p);
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p); 

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

//  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  Spin spin(sp); 
  mode::Ptr mode(mode::make(pplm ,qmp));

  JaynesCummings<> jc(spin,mode,pjc);

  // BinarySystem system(jc);

  StateVector1 psiMode(mode::init(pplm)), psiSpin(spin.getDimension());

  psiSpin()(0)=1; psiSpin.renorm();

  StateVector2 psi(psiSpin*psiMode);

//  StateVector psi(qbit::init(pplqb)*mode::init(pplm));
  // entangled state: 
  // psi(qbit::state0()*mode::fock(1,pplm.cutoff)+qbit::state1()*mode::fock(0,pplm.cutoff));
  /*()+qbit::state1()())*
  psi()(0,0)=1; psi()(1,0)=1;
  */
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,binary::make(jc),pe);



  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
