// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"
#include "Mode.h"

#include "SmartPtr.h"

using cpputils::nonOwningConstSharedPtr;


struct DummyQuaternary : structure::Interaction<4>
{
  DummyQuaternary(mode::Ptr m0, mode::Ptr m1, mode::Ptr m2, mode::Ptr m3) 
    : structure::Interaction<4>(Frees(m0,m1,m2,m3)) {}
};

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars pe(p);

  mode::ParsPumpedLossy
    pm0(p,"0"),
    pm1(p,"1"),
    pm2(p,"2"),
    pm3(p,"3");

  update(p,argc,argv,"--");

  PumpedLossyMode<> 
    m0(pm0), m2(pm2);

  PumpedLossyModeAlternative<false>
    m1(pm1), m3(pm3);

  DummyQuaternary dq(nonOwningConstSharedPtr(&m0),nonOwningConstSharedPtr(&m1),nonOwningConstSharedPtr(&m2),nonOwningConstSharedPtr(&m3));

  quantumdata::StateVector<4> psi(init(pm0)*init(pm1)*init(pm2)*init(pm3));


  evolve(psi,
         composite::make(
                         _<0,1,2,3>(dq)
                        ),
         pe);


}
