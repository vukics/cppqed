#include "EvolutionBinary.h"

#include "JaynesCummings.h"

#include "MatrixOfHamiltonian.h"


using namespace std;

typedef quantumdata::StateVector<2> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  evolution::Pars pe(p); // Driver Parameters
  qbit::ParsPumpedLossy pplqb(p); 
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p); 

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  mode::Ptr mode(mode::make(pplm ,qmp/*,structure::averaged::ReducedDensityOperator("Mode",pplm.cutoff,true)*/));

  // BinarySystem system(jc);
  /*
  mode::StateVector psiMode(pplm.cutoff);
  psiMode=(mode::fock(0,pplm.cutoff,mathutils::PI/2.)+mode::fock(1,pplm.cutoff)+mode::fock(2,pplm.cutoff,-mathutils::PI/2.))/sqrt(3.);
  */
  StateVector psi(qbit::init(pplqb)*mode::init(pplm));
  // entangled state: 
  // psi=qbit::state0()*mode::fock(1,pplm.cutoff)+qbit::state1()*mode::fock(0,pplm.cutoff);
  psi.renorm();

  evolve(psi,binary::make(jaynescummings::make(qbit,mode,pjc)),pe);






}
