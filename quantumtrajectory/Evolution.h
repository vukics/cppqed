// -*- C++ -*-
#ifndef _EVOLUTION_MODE_H
#define _EVOLUTION_MODE_H

#include "EvolutionFwd.h"

#include "StateVectorFwd.h"
#include "QuantumSystemFwd.h"

#include "TMP_Tools.h"

#include<iosfwd>

std::ostream& operator<<(std::ostream&, EvolutionMode );
std::istream& operator>>(std::istream&, EvolutionMode&);


template<int RANK, typename V>
void evolve(quantumdata::StateVector<RANK>&, const structure::QuantumSystem<RANK>&,
	    const ParsEvolution&,
	    V);
// NEEDS_WORK a version for nonorthogonals should also be provided...


template<int RANK>
inline
void evolve(quantumdata::StateVector<RANK>& psi,
	    const structure::QuantumSystem<RANK>& sys,
	    const ParsEvolution& p)
{
  evolve(psi,sys,p,tmptools::V0());
}


#include "impl/Evolution.tcc"

#endif // _EVOLUTION_MODE_H
