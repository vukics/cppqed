#include "ParsEvolution.h"
#include "Simulated.h"

#include "PumpedTwoLevelAtom.h"
#include "Qbit.h"

#include "impl/StateVector.tcc"
#include "impl/DensityOperator.tcc"

using namespace std       ;
using namespace cpputils  ;
using namespace trajectory;
using namespace qbit      ;

typedef TTD_DARRAY(1) Array;

void derivs(double, const Array& b, Array& dbdt, const ParsPumpedLossy& p)
{
  double
    z=b(0);

  dcomp
    Omega(-p.gamma,p.delta),
    s(b(1),-b(2));

  dcomp
    temp(-2.*conj(p.eta)*z+conj(Omega)*s);

  dbdt(0)= 2.*real(p.eta*s)+2*p.gamma*(1-z);
  dbdt(1)= real(temp);
  dbdt(2)=-imag(temp);
}


int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsEvolution pt(p);
  ParsPumpedLossy pp2la(p); 

  // Parameter finalization
  update(p,argc,argv,"--");

  PumpedTwoLevelAtomSch atom(pp2la);

  double dtinit=.1/atom.highestFrequency();
 
  Array zxy(3);

  {  
    quantumdata::DensityOperator<1> rho(qbit::init(pp2la));

    zxy(0)= 2*real(rho()(0,0))-1;
    zxy(1)= 2*real(rho()(0,1));
    zxy(2)=-2*imag(rho()(0,1));
  }

  Simulated<Array> S(zxy,bind(derivs,_1,_2,_3,pp2la),dtinit,Array(),pt);
  
  evolve(S,pt);
  
}
