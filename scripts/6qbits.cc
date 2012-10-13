#define BLITZ_ARRAY_LARGEST_RANK 15

#include "EvolutionComposite.h"
#include "Mode.h"

#include "JaynesCummings.h"

//#include "BinarySystem.h"

#include "MatrixOfHamiltonian.h"


using namespace std;

typedef quantumdata::StateVector<7> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  qbit::ParsPumpedLossy pplqb(p);
  qbit::ParsPumpedLossy pplqb2(p,"2");
  qbit::ParsPumpedLossy pplqb3(p,"3"); 
  qbit::ParsPumpedLossy pplqb4(p,"4");
  qbit::ParsPumpedLossy pplqb5(p,"5");
  qbit::ParsPumpedLossy pplqb6(p,"6"); 
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);
  jaynescummings::Pars  pjc2  (p,"2");
  jaynescummings::Pars  pjc3  (p,"3");
  jaynescummings::Pars  pjc4  (p,"4"); 
  jaynescummings::Pars  pjc5  (p,"5");
  jaynescummings::Pars  pjc6  (p,"6");

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  qbit::Ptr qbit2(qbit::make(pplqb2,qmp));
  qbit::Ptr qbit3(qbit::make(pplqb3,qmp));
  qbit::Ptr qbit4(qbit::make(pplqb4,qmp));
  qbit::Ptr qbit5(qbit::make(pplqb5,qmp));
  qbit::Ptr qbit6(qbit::make(pplqb6,qmp));

  mode::Ptr mode(mode::make(pplm ,qmp));

  JaynesCummings<> jc(qbit,mode,pjc);
  JaynesCummings<> jc2(qbit2,mode,pjc2);
  JaynesCummings<> jc3(qbit3,mode,pjc3);
  JaynesCummings<> jc4(qbit4,mode,pjc4);
  JaynesCummings<> jc5(qbit5,mode,pjc5);
  JaynesCummings<> jc6(qbit6,mode,pjc6);

  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb2)*qbit::init(pplqb3)*qbit::init(pplqb4)*qbit::init(pplqb5)*qbit::init(pplqb6)*mode::init(pplm));
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,composite::make(Act<0,6>(jc),Act<1,6>(jc2),Act<2,6>(jc3),Act<3,6>(jc4),Act<4,6>(jc5),Act<5,6>(jc6)),pe);
 
  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
