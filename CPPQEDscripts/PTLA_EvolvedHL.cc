// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution_.h"

#include "StateVector.h"
#include "DensityOperator.h"
#include "PumpedTwoLevelAtom.h"
#include "Qbit.h"

#include "Simulated.h"

#include<algorithm>

using namespace std       ;
using namespace cppqedutils  ;
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
 
  Array sigma(2);

  {
    quantumdata::DensityOperator<1> rho(qbit::init(pp2la));

    sigma=real(rho(1)(0)),imag(rho(1)(0));
  }

  Simulated<Array> S(sigma,[&](double, const Array& b, Array& dbdt) {
    dcomp
      Omega(-pp2la.gamma,pp2la.delta),
      sigma(b(0),b(1));

    dcomp
      temp(Omega*sigma-pp2la.eta);

    dbdt=real(temp),imag(temp);
  },dtinit,pt);
  
  run(S,pt);

  
}
