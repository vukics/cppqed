// -*- C++ -*-
#ifndef _EVOLUTION_MODE_H
#define _EVOLUTION_MODE_H

#include "EvolutionFwd.h"

#include "StateVectorFwd.h"
#include "QuantumSystemFwd.h"

#include "TMP_Tools.h"

#include <boost/shared_ptr.hpp>

#include <iosfwd>


std::ostream& operator<<(std::ostream&, EvolutionMode );
std::istream& operator>>(std::istream&, EvolutionMode&);


template<int RANK, typename V>
void evolve(quantumdata::StateVector<RANK>&, const structure::QuantumSystem<RANK>&,
	    const ParsEvolution&,
	    V);



template<int RANK>
inline
void evolve(quantumdata::StateVector<RANK>& psi,
	    const structure::QuantumSystem<RANK>& sys,
	    const ParsEvolution& p)
{
  evolve(psi,sys,p,tmptools::V_Empty());
}


template<int RANK, typename SYS, typename V>
inline
void evolve(quantumdata::StateVector<RANK>& psi,
	    boost::shared_ptr<const SYS> sys,
	    const ParsEvolution& p,
	    V v)
{
  evolve(psi,*sys,p,v);
}


template<int RANK, typename SYS>
inline
void evolve(quantumdata::StateVector<RANK>& psi,
	    boost::shared_ptr<const SYS> sys,
	    const ParsEvolution& p)
{
  evolve(psi,*sys,p);
}


// C++11: use default template argument tmptools::V_Empty to fuse the last two functions


#include "impl/Evolution.tcc"

#endif // _EVOLUTION_MODE_H
