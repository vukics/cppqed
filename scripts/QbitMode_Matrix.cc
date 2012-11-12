#include "EvolutionBinary.h"
#include "Mode.h"

#include "TimeIndependentMatrixHamiltonian.h"
#include "JaynesCummings.h"

#include "MatrixOfHamiltonian.h"


using namespace std;
using namespace blitz;

using linalg::CMatrix;

typedef quantumdata::StateVector<1> StateVector1;
typedef StateVector1::StateVectorLow StateVectorLow1;

typedef quantumdata::StateVector<2> StateVector2;
typedef StateVector2::StateVectorLow StateVectorLow2;

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
  
  pe.noise=false;

  qmp=QMP_SCH; // Important!!! MatrixHamiltonian only works in this case.

  // ****** ****** ****** ****** ****** ******

  binary::Ptr system2(binary::make(jaynescummings::make(qbit::make(pplqb,qmp),mode::make(pplm ,qmp),pjc)));

  StateVector2 psi2(qbit::init(pplqb)*mode::init(pplm)); psi2.renorm();
  StateVector1 psi1(StateVectorLow1(psi2().data(),shape(system2->getTotalDimension()),neverDeleteData),quantumdata::byReference);
  // Now psi1 and psi2 referece the same data.

  const CMatrix hamiltonian(calculateMatrix(*system2));

  evolve(psi1,TimeIndependentMatrixHamiltonianAveraged<2,true>(CMatrix(hamiltonian/DCOMP_I),system2,psi2),pe);
  
  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
