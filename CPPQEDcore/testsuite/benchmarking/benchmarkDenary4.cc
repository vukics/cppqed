// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "benchmarking.h"

#include "Mode_.h"

using namespace mode;


struct Ia : structure::Interaction<4>, structure::TridiagonalHamiltonian<4,false>
{
  Ia(const ModeBase* m) : structure::Interaction<4>(Frees(m,m,m,m)), structure::TridiagonalHamiltonian<4,false>(aop(m)*aop(m).dagger()*aop(m)*aop(m).dagger()) {}
};


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  const size_t cutoff=4;

  // ****** ****** ****** ****** ****** ******

  const ModeBase mode(cutoff);

  const Ia ia(&mode);
  const structure::Interaction<10> dummy(structure::Interaction<10>::Frees(&mode,&mode,&mode,&mode,&mode,&mode,&mode,&mode,&mode,&mode));

  const Composite<composite::result_of::make_vector<Act<0,2,7,8>,Act<0,1,2,3,4,5,6,7,8,9> >::type>
    sys(makeComposite(Act<0,2,7,8>(ia),
                      Act<0,1,2,3,4,5,6,7,8,9>(dummy)
                                     ));

  benchmark(sys,ia,Vector<0,2,7,8>());

}




namespace quantumoperator {


template<> template<>
void Tridiagonal<4>::apply<4>(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  using blitz::Range;
  
  const Diagonals & a=diagonals_  ;
  const Dimensions& K=differences_;

  Range
    I0 (0,psi.ubound(0)),
    I0l(0,psi.ubound(0)-int(K(0))),
    I0h(I0l+int(K(0))),
    I1 (0,psi.ubound(1)),
    I1l(0,psi.ubound(1)-int(K(1))),
    I1h(I1l+int(K(1))),
    I2 (0,psi.ubound(2)),
    I2l(0,psi.ubound(2)-int(K(2))),
    I2h(I2l+int(K(2))),
    I3 (0,psi.ubound(3)),
    I3l(0,psi.ubound(3)-int(K(3))),
    I3h(I3l+int(K(3)));

  if (a(0).size()) dpsidt(I0 ,I1 ,I2 ,I3 )+=a(0)(I0 ,I1 ,I2 ,I3 )*psi(I0 ,I1 ,I2 ,I3 );
  if (a(1).size()) dpsidt(I0 ,I1 ,I2 ,I3h)+=a(1)(I0 ,I1 ,I2 ,I3l)*psi(I0 ,I1 ,I2 ,I3l);
  if (a(2).size()) dpsidt(I0 ,I1 ,I2 ,I3l)+=a(2)(I0 ,I1 ,I2 ,I3l)*psi(I0 ,I1 ,I2 ,I3h);
  if (a(3).size()) dpsidt(I0 ,I1 ,I2h,I3 )+=a(3)(I0 ,I1 ,I2l,I3 )*psi(I0 ,I1 ,I2l,I3 );
  if (a(4).size()) dpsidt(I0 ,I1 ,I2h,I3h)+=a(4)(I0 ,I1 ,I2l,I3l)*psi(I0 ,I1 ,I2l,I3l);
  if (a(5).size()) dpsidt(I0 ,I1 ,I2h,I3l)+=a(5)(I0 ,I1 ,I2l,I3l)*psi(I0 ,I1 ,I2l,I3h);
  if (a(6).size()) dpsidt(I0 ,I1 ,I2l,I3 )+=a(6)(I0 ,I1 ,I2l,I3 )*psi(I0 ,I1 ,I2h,I3 );
  if (a(7).size()) dpsidt(I0 ,I1 ,I2l,I3h)+=a(7)(I0 ,I1 ,I2l,I3l)*psi(I0 ,I1 ,I2h,I3l);
  if (a(8).size()) dpsidt(I0 ,I1 ,I2l,I3l)+=a(8)(I0 ,I1 ,I2l,I3l)*psi(I0 ,I1 ,I2h,I3h);
  if (a(9).size()) dpsidt(I0 ,I1h,I2 ,I3 )+=a(9)(I0 ,I1l,I2 ,I3 )*psi(I0 ,I1l,I2 ,I3 );
  if (a(10).size()) dpsidt(I0 ,I1h,I2 ,I3h)+=a(10)(I0 ,I1l,I2 ,I3l)*psi(I0 ,I1l,I2 ,I3l);
  if (a(11).size()) dpsidt(I0 ,I1h,I2 ,I3l)+=a(11)(I0 ,I1l,I2 ,I3l)*psi(I0 ,I1l,I2 ,I3h);
  if (a(12).size()) dpsidt(I0 ,I1h,I2h,I3 )+=a(12)(I0 ,I1l,I2l,I3 )*psi(I0 ,I1l,I2l,I3 );
  if (a(13).size()) dpsidt(I0 ,I1h,I2h,I3h)+=a(13)(I0 ,I1l,I2l,I3l)*psi(I0 ,I1l,I2l,I3l);
  if (a(14).size()) dpsidt(I0 ,I1h,I2h,I3l)+=a(14)(I0 ,I1l,I2l,I3l)*psi(I0 ,I1l,I2l,I3h);
  if (a(15).size()) dpsidt(I0 ,I1h,I2l,I3 )+=a(15)(I0 ,I1l,I2l,I3 )*psi(I0 ,I1l,I2h,I3 );
  if (a(16).size()) dpsidt(I0 ,I1h,I2l,I3h)+=a(16)(I0 ,I1l,I2l,I3l)*psi(I0 ,I1l,I2h,I3l);
  if (a(17).size()) dpsidt(I0 ,I1h,I2l,I3l)+=a(17)(I0 ,I1l,I2l,I3l)*psi(I0 ,I1l,I2h,I3h);
  if (a(18).size()) dpsidt(I0 ,I1l,I2 ,I3 )+=a(18)(I0 ,I1l,I2 ,I3 )*psi(I0 ,I1h,I2 ,I3 );
  if (a(19).size()) dpsidt(I0 ,I1l,I2 ,I3h)+=a(19)(I0 ,I1l,I2 ,I3l)*psi(I0 ,I1h,I2 ,I3l);
  if (a(20).size()) dpsidt(I0 ,I1l,I2 ,I3l)+=a(20)(I0 ,I1l,I2 ,I3l)*psi(I0 ,I1h,I2 ,I3h);
  if (a(21).size()) dpsidt(I0 ,I1l,I2h,I3 )+=a(21)(I0 ,I1l,I2l,I3 )*psi(I0 ,I1h,I2l,I3 );
  if (a(22).size()) dpsidt(I0 ,I1l,I2h,I3h)+=a(22)(I0 ,I1l,I2l,I3l)*psi(I0 ,I1h,I2l,I3l);
  if (a(23).size()) dpsidt(I0 ,I1l,I2h,I3l)+=a(23)(I0 ,I1l,I2l,I3l)*psi(I0 ,I1h,I2l,I3h);
  if (a(24).size()) dpsidt(I0 ,I1l,I2l,I3 )+=a(24)(I0 ,I1l,I2l,I3 )*psi(I0 ,I1h,I2h,I3 );
  if (a(25).size()) dpsidt(I0 ,I1l,I2l,I3h)+=a(25)(I0 ,I1l,I2l,I3l)*psi(I0 ,I1h,I2h,I3l);
  if (a(26).size()) dpsidt(I0 ,I1l,I2l,I3l)+=a(26)(I0 ,I1l,I2l,I3l)*psi(I0 ,I1h,I2h,I3h);
  if (a(27).size()) dpsidt(I0h,I1 ,I2 ,I3 )+=a(27)(I0l,I1 ,I2 ,I3 )*psi(I0l,I1 ,I2 ,I3 );
  if (a(28).size()) dpsidt(I0h,I1 ,I2 ,I3h)+=a(28)(I0l,I1 ,I2 ,I3l)*psi(I0l,I1 ,I2 ,I3l);
  if (a(29).size()) dpsidt(I0h,I1 ,I2 ,I3l)+=a(29)(I0l,I1 ,I2 ,I3l)*psi(I0l,I1 ,I2 ,I3h);
  if (a(30).size()) dpsidt(I0h,I1 ,I2h,I3 )+=a(30)(I0l,I1 ,I2l,I3 )*psi(I0l,I1 ,I2l,I3 );
  if (a(31).size()) dpsidt(I0h,I1 ,I2h,I3h)+=a(31)(I0l,I1 ,I2l,I3l)*psi(I0l,I1 ,I2l,I3l);
  if (a(32).size()) dpsidt(I0h,I1 ,I2h,I3l)+=a(32)(I0l,I1 ,I2l,I3l)*psi(I0l,I1 ,I2l,I3h);
  if (a(33).size()) dpsidt(I0h,I1 ,I2l,I3 )+=a(33)(I0l,I1 ,I2l,I3 )*psi(I0l,I1 ,I2h,I3 );
  if (a(34).size()) dpsidt(I0h,I1 ,I2l,I3h)+=a(34)(I0l,I1 ,I2l,I3l)*psi(I0l,I1 ,I2h,I3l);
  if (a(35).size()) dpsidt(I0h,I1 ,I2l,I3l)+=a(35)(I0l,I1 ,I2l,I3l)*psi(I0l,I1 ,I2h,I3h);
  if (a(36).size()) dpsidt(I0h,I1h,I2 ,I3 )+=a(36)(I0l,I1l,I2 ,I3 )*psi(I0l,I1l,I2 ,I3 );
  if (a(37).size()) dpsidt(I0h,I1h,I2 ,I3h)+=a(37)(I0l,I1l,I2 ,I3l)*psi(I0l,I1l,I2 ,I3l);
  if (a(38).size()) dpsidt(I0h,I1h,I2 ,I3l)+=a(38)(I0l,I1l,I2 ,I3l)*psi(I0l,I1l,I2 ,I3h);
  if (a(39).size()) dpsidt(I0h,I1h,I2h,I3 )+=a(39)(I0l,I1l,I2l,I3 )*psi(I0l,I1l,I2l,I3 );
  if (a(40).size()) dpsidt(I0h,I1h,I2h,I3h)+=a(40)(I0l,I1l,I2l,I3l)*psi(I0l,I1l,I2l,I3l);
  if (a(41).size()) dpsidt(I0h,I1h,I2h,I3l)+=a(41)(I0l,I1l,I2l,I3l)*psi(I0l,I1l,I2l,I3h);
  if (a(42).size()) dpsidt(I0h,I1h,I2l,I3 )+=a(42)(I0l,I1l,I2l,I3 )*psi(I0l,I1l,I2h,I3 );
  if (a(43).size()) dpsidt(I0h,I1h,I2l,I3h)+=a(43)(I0l,I1l,I2l,I3l)*psi(I0l,I1l,I2h,I3l);
  if (a(44).size()) dpsidt(I0h,I1h,I2l,I3l)+=a(44)(I0l,I1l,I2l,I3l)*psi(I0l,I1l,I2h,I3h);
  if (a(45).size()) dpsidt(I0h,I1l,I2 ,I3 )+=a(45)(I0l,I1l,I2 ,I3 )*psi(I0l,I1h,I2 ,I3 );
  if (a(46).size()) dpsidt(I0h,I1l,I2 ,I3h)+=a(46)(I0l,I1l,I2 ,I3l)*psi(I0l,I1h,I2 ,I3l);
  if (a(47).size()) dpsidt(I0h,I1l,I2 ,I3l)+=a(47)(I0l,I1l,I2 ,I3l)*psi(I0l,I1h,I2 ,I3h);
  if (a(48).size()) dpsidt(I0h,I1l,I2h,I3 )+=a(48)(I0l,I1l,I2l,I3 )*psi(I0l,I1h,I2l,I3 );
  if (a(49).size()) dpsidt(I0h,I1l,I2h,I3h)+=a(49)(I0l,I1l,I2l,I3l)*psi(I0l,I1h,I2l,I3l);
  if (a(50).size()) dpsidt(I0h,I1l,I2h,I3l)+=a(50)(I0l,I1l,I2l,I3l)*psi(I0l,I1h,I2l,I3h);
  if (a(51).size()) dpsidt(I0h,I1l,I2l,I3 )+=a(51)(I0l,I1l,I2l,I3 )*psi(I0l,I1h,I2h,I3 );
  if (a(52).size()) dpsidt(I0h,I1l,I2l,I3h)+=a(52)(I0l,I1l,I2l,I3l)*psi(I0l,I1h,I2h,I3l);
  if (a(53).size()) dpsidt(I0h,I1l,I2l,I3l)+=a(53)(I0l,I1l,I2l,I3l)*psi(I0l,I1h,I2h,I3h);
  if (a(54).size()) dpsidt(I0l,I1 ,I2 ,I3 )+=a(54)(I0l,I1 ,I2 ,I3 )*psi(I0h,I1 ,I2 ,I3 );
  if (a(55).size()) dpsidt(I0l,I1 ,I2 ,I3h)+=a(55)(I0l,I1 ,I2 ,I3l)*psi(I0h,I1 ,I2 ,I3l);
  if (a(56).size()) dpsidt(I0l,I1 ,I2 ,I3l)+=a(56)(I0l,I1 ,I2 ,I3l)*psi(I0h,I1 ,I2 ,I3h);
  if (a(57).size()) dpsidt(I0l,I1 ,I2h,I3 )+=a(57)(I0l,I1 ,I2l,I3 )*psi(I0h,I1 ,I2l,I3 );
  if (a(58).size()) dpsidt(I0l,I1 ,I2h,I3h)+=a(58)(I0l,I1 ,I2l,I3l)*psi(I0h,I1 ,I2l,I3l);
  if (a(59).size()) dpsidt(I0l,I1 ,I2h,I3l)+=a(59)(I0l,I1 ,I2l,I3l)*psi(I0h,I1 ,I2l,I3h);
  if (a(60).size()) dpsidt(I0l,I1 ,I2l,I3 )+=a(60)(I0l,I1 ,I2l,I3 )*psi(I0h,I1 ,I2h,I3 );
  if (a(61).size()) dpsidt(I0l,I1 ,I2l,I3h)+=a(61)(I0l,I1 ,I2l,I3l)*psi(I0h,I1 ,I2h,I3l);
  if (a(62).size()) dpsidt(I0l,I1 ,I2l,I3l)+=a(62)(I0l,I1 ,I2l,I3l)*psi(I0h,I1 ,I2h,I3h);
  if (a(63).size()) dpsidt(I0l,I1h,I2 ,I3 )+=a(63)(I0l,I1l,I2 ,I3 )*psi(I0h,I1l,I2 ,I3 );
  if (a(64).size()) dpsidt(I0l,I1h,I2 ,I3h)+=a(64)(I0l,I1l,I2 ,I3l)*psi(I0h,I1l,I2 ,I3l);
  if (a(65).size()) dpsidt(I0l,I1h,I2 ,I3l)+=a(65)(I0l,I1l,I2 ,I3l)*psi(I0h,I1l,I2 ,I3h);
  if (a(66).size()) dpsidt(I0l,I1h,I2h,I3 )+=a(66)(I0l,I1l,I2l,I3 )*psi(I0h,I1l,I2l,I3 );
  if (a(67).size()) dpsidt(I0l,I1h,I2h,I3h)+=a(67)(I0l,I1l,I2l,I3l)*psi(I0h,I1l,I2l,I3l);
  if (a(68).size()) dpsidt(I0l,I1h,I2h,I3l)+=a(68)(I0l,I1l,I2l,I3l)*psi(I0h,I1l,I2l,I3h);
  if (a(69).size()) dpsidt(I0l,I1h,I2l,I3 )+=a(69)(I0l,I1l,I2l,I3 )*psi(I0h,I1l,I2h,I3 );
  if (a(70).size()) dpsidt(I0l,I1h,I2l,I3h)+=a(70)(I0l,I1l,I2l,I3l)*psi(I0h,I1l,I2h,I3l);
  if (a(71).size()) dpsidt(I0l,I1h,I2l,I3l)+=a(71)(I0l,I1l,I2l,I3l)*psi(I0h,I1l,I2h,I3h);
  if (a(72).size()) dpsidt(I0l,I1l,I2 ,I3 )+=a(72)(I0l,I1l,I2 ,I3 )*psi(I0h,I1h,I2 ,I3 );
  if (a(73).size()) dpsidt(I0l,I1l,I2 ,I3h)+=a(73)(I0l,I1l,I2 ,I3l)*psi(I0h,I1h,I2 ,I3l);
  if (a(74).size()) dpsidt(I0l,I1l,I2 ,I3l)+=a(74)(I0l,I1l,I2 ,I3l)*psi(I0h,I1h,I2 ,I3h);
  if (a(75).size()) dpsidt(I0l,I1l,I2h,I3 )+=a(75)(I0l,I1l,I2l,I3 )*psi(I0h,I1h,I2l,I3 );
  if (a(76).size()) dpsidt(I0l,I1l,I2h,I3h)+=a(76)(I0l,I1l,I2l,I3l)*psi(I0h,I1h,I2l,I3l);
  if (a(77).size()) dpsidt(I0l,I1l,I2h,I3l)+=a(77)(I0l,I1l,I2l,I3l)*psi(I0h,I1h,I2l,I3h);
  if (a(78).size()) dpsidt(I0l,I1l,I2l,I3 )+=a(78)(I0l,I1l,I2l,I3 )*psi(I0h,I1h,I2h,I3 );
  if (a(79).size()) dpsidt(I0l,I1l,I2l,I3h)+=a(79)(I0l,I1l,I2l,I3l)*psi(I0h,I1h,I2h,I3l);
  if (a(80).size()) dpsidt(I0l,I1l,I2l,I3l)+=a(80)(I0l,I1l,I2l,I3l)*psi(I0h,I1h,I2h,I3h);

}


} // quantumoperator

