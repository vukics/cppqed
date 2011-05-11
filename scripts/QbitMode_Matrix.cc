#include "EvolutionHigh.h"

#include "TimeIndependentMatrixHamiltonian.h"
#include "JaynesCummings.h"

#include "BinarySystem.h"

#include "MatrixOfHamiltonian.h"


using namespace std;
using namespace blitz;

using formdouble::special;
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

  qbit::SmartPtr qbit(qbit::maker(pplqb,qmp));
  mode::SmartPtr mode(mode::maker(pplm ,qmp));

  JaynesCummings<> jc(qbit,mode,pjc);

  BinarySystem system2(jc);

  StateVector2 psi2(qbit::init(pplqb)*mode::init(pplm)); psi2.renorm();
  StateVector1 psi1(StateVectorLow1(psi2().data(),shape(system2.getTotalDimension()),neverDeleteData),quantumdata::byReference);
  // Now psi1 and psi2 referece the same data.

  const CMatrix hamiltonian(calculateMatrix(system2));

  TimeIndependentMatrixHamiltonian system1(CMatrix(hamiltonian/DCOMP_I));

  quantumtrajectory::MCWF_Trajectory<1> traj(psi1,system1,pe);
  
  {

#define TRAJ_DISPLAY traj.getOstream()<<special()(traj.getTime())<<special()(traj.getDtDid()); structure::Averaged<2>::display(0,psi2,traj.getOstream(),pe.precision,&system2); traj.getOstream()<<endl;

    traj.displayParameters();
    traj.getOstream()<<"# Run Trajectory. Displaying in every "<<pe.Dt<<endl<<endl
		     <<"# Key to data:"<<endl;
    traj.displayKey();
    traj.getOstream()<<endl;
    TRAJ_DISPLAY
    while (traj.getTime()<pe.T) {
      traj.evolve(std::min(pe.Dt,pe.T-traj.getTime()));
      TRAJ_DISPLAY
    }

  }


  } catch (ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
