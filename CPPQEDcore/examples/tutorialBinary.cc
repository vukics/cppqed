#include "EvolutionBinary.h"

#include "JaynesCummings.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsEvolution pe(p); 

  qbit::ParsPumpedLossy pq(p);
  mode::ParsPumpedLossy pm(p);
  jaynescummings::Pars pjc(p);

  update(p,argc,argv,"--"); // Parsing the command line

  quantumdata::StateVector<2> psi(init(pq)*init(pm));
  evolve(psi,
         binary::make(JaynesCummings<>(make(pq,QMP_IP),
                                       make(pm,QMP_IP),
                                       pjc)),
         pe);

}
