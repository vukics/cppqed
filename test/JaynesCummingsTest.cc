#include "JaynesCummings.h"

#include "Mode.h"
#include "Qbit.h"

#include "StateVector.h"

#include "Pars.h"

using namespace std;

typedef quantumdata::StateVector<2> StateVector;
typedef StateVector::Dimensions     Dimensions ;

int main()
{
  parameters::ParameterTable p;
  ParsJaynesCummings  pjc  (p); 
  pjc.g=dcomp(1,2);

  ModeBase mb(5,0);
  QbitBase qb;
  JaynesCummings<false> jc(qb,mb,pjc);

  StateVector psi(Dimensions(2,5));

  psi()(1,0)=1;
  cout<<psi();

  StateVector dpsidt(psi.getDimensions());
  structure::Hamiltonian<2>::addContribution(0,psi(),dpsidt(),0,&jc);
  cout<<dpsidt();

}
