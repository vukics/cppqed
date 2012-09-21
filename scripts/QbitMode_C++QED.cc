#include "Evolution.h"
#include "Mode.h"

#include "JaynesCummings.h"

#include "BinarySystem.h"

#include "MatrixOfHamiltonian.h"


using namespace std;

typedef quantumdata::StateVector<2> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  qbit::ParsPumpedLossy pplqb(p); 
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p); 

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  mode::Ptr mode(mode::make(pplm ,qmp));

  JaynesCummings<> jc(qbit,mode,pjc);

  // BinarySystem system(jc);

  StateVector psi(qbit::init(pplqb)*mode::init(pplm));
  // entangled state: 
  // psi(qbit::state0()*mode::fock(1,pplm.cutoff)+qbit::state1()*mode::fock(0,pplm.cutoff));
  /*()+qbit::state1()())*
  psi()(0,0)=1; psi()(1,0)=1;
  */
  psi.renorm();

  evolve(psi,binary::make(jc),pe,tmptools::Vector<0>());


  // The 3 further ways to create a JaynesCummings:
  {
    PumpedLossyQbitSch qbitv(pplqb);
    PumpedLossyMode<>  modev(pplm );
    {
      JaynesCummings<> jc(qbitv,mode ,pjc);
    }
    {
      JaynesCummings<> jc(qbit ,modev,pjc);
    }
    {
      JaynesCummings<> jc(qbitv,modev,pjc);
    }
  }



  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
