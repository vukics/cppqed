// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution_.h"

#include "PumpedTwoLevelAtom.h"
#include "Qbit.h"

#include "StateVector.h"
#include "DensityOperator.h"

#include "Simulated.h"

using namespace std       ;
using namespace cpputils  ;
using namespace trajectory;
using namespace qbit      ;

typedef DArray<1> Array;


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pt(p);
  ParsPumpedLossy pp2la(p); 

  // Parameter finalization
  update(p,argc,argv,"--");

  PumpedTwoLevelAtomSch atom(pp2la);

  double dtinit=.1/atom.highestFrequency();
 
  Array zxy(3);

  {  
    quantumdata::DensityOperator<1> rho(qbit::init(pp2la));

    zxy=
      2*real(rho(0)(0))-1,
      2*real(rho(0)(1))  ,
     -2*imag(rho(0)(1))  ;
  }

  Simulated<Array> S(zxy,[&](double, const Array& b, Array& dbdt) {
    double z=b(0);

    dcomp Omega(-pp2la.gamma,pp2la.delta), s(b(1),-b(2));

    dcomp temp(-2.*conj(pp2la.eta)*z+conj(Omega)*s);

    dbdt=
      2.*real(pp2la.eta*s)+2*pp2la.gamma*(1-z),
       real(temp),
      -imag(temp);    
  },dtinit,pt);
  
  run(S,pt);
  
}
