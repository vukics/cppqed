// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution_.h"
#include "Simulated.h"

#include "StateVector.tcc"
#include "DensityOperator.tcc"
#include "PumpedTwoLevelAtom.h"
#include "Qbit.h"

#include<algorithm>

using namespace std       ;
using namespace cpputils  ;
using namespace trajectory;
using namespace qbit      ;

typedef DArray<1> Array;

void derivs(double, const Array& b, Array& dbdt, const ParsPumpedLossy& p)
{
  dcomp
    Omega(-p.gamma,p.delta),
    sigma(b(0),b(1));

  dcomp
    temp(Omega*sigma-p.eta);

  dbdt=real(temp),imag(temp);
}


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars pt(p);
  ParsPumpedLossy pp2la(p); 

  // Parameter finalization
  update(p,argc,argv,"--");

  PumpedTwoLevelAtomSch atom(pp2la);

  double dtinit=.1/atom.highestFrequency();
 
  Array sigma(2);

  {
    quantumdata::DensityOperator<1> rho(qbit::init(pp2la));

    sigma=real(rho(1)(0)),imag(rho(1)(0));
  }

  Simulated<Array> S(sigma,bind(derivs,_1,_2,_3,pp2la),dtinit,pt);
  
  run(S,pt);

  
}
