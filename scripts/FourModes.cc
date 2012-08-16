#include "Evolution.h"
#include "Mode.h"

#include "Composite.h"


struct DummyQuaternary : structure::Interaction<4>
{
  DummyQuaternary(ModeBase* m0, ModeBase* m1, ModeBase* m2, ModeBase* m3) 
    : structure::Interaction<4>(Frees(m0,m1,m2,m3)) {}
};

int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsEvolution pe(p);

  mode::ParsPumpedLossy pm0(p,"0");
  mode::ParsPumpedLossy pm1(p,"1");
  mode::ParsPumpedLossy pm2(p,"2");
  mode::ParsPumpedLossy pm3(p,"3");

  update(p,argc,argv,"--");

  PumpedLossyMode<> 
    m0(pm0), m2(pm2);

  PumpedLossyModeAlternative<>
    m1(pm1), m3(pm3);

  DummyQuaternary dq(&m0,&m1,&m2,&m3);

  quantumdata::StateVector<4> psi(init(pm0)*init(pm1)*init(pm2)*init(pm3));


  evolve(psi,
	 composite::make(	         
			 Act<0,1,2,3>(dq)
					 ),
	 pe);


}
