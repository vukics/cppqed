// -*- C++ -*-
#ifndef STRUCTURE_STRUCTURE_H_INCLUDED
#define STRUCTURE_STRUCTURE_H_INCLUDED

#include "QuantumSystemFwd.h"
#include "DynamicsBaseFwd.h"

#include "Exact.h"
#include "Hamiltonian.h"
#include "Liouvillean.h"
#include "Averaged.h"


// Note that, in each case, the most general class is taken.
// This is because Hamiltonian ONE_TIME & NO_TIME is derived from TWO_TIME, Liouvillean false is derived from true, and Averaged false is derived from true.
// At the same time, these most general cases are the default ones.

namespace structure {

using boost::dynamic_pointer_cast;

template<int RANK>
inline 
const typename Exact<RANK>::Ptr 
qse(typename QuantumSystem<RANK>::Ptr quantumSystem)
{return dynamic_pointer_cast<const Exact<RANK> >(quantumSystem);}

template<int RANK>
inline 
const typename Hamiltonian<RANK,TWO_TIME>::Ptr 
qsh(typename QuantumSystem<RANK>::Ptr quantumSystem)
{return dynamic_pointer_cast<const Hamiltonian<RANK,TWO_TIME> >(quantumSystem);}

template<int RANK>
inline 
const typename Liouvillean<RANK,true>::Ptr 
qsl(typename QuantumSystem<RANK>::Ptr quantumSystem)
{return dynamic_pointer_cast<const Liouvillean<RANK,true> >(quantumSystem);}

template<int RANK>
inline 
const typename Averaged<RANK,true>::Ptr 
qsa(typename QuantumSystem<RANK>::Ptr quantumSystem)
{return dynamic_pointer_cast<const Averaged<RANK,true> >(quantumSystem);}



template<int RANK>
inline 
const typename Exact<RANK>::Ptr 
qse(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Exact<RANK> >(base);}

template<int RANK>
inline 
const typename Hamiltonian<RANK,TWO_TIME>::Ptr 
qsh(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Hamiltonian<RANK,TWO_TIME> >(base);}

template<int RANK>
inline 
const typename Liouvillean<RANK,true>::Ptr 
qsl(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Liouvillean<RANK,true> >(base);}

template<int RANK>
inline 
const typename Averaged<RANK,true>::Ptr 
qsa(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Averaged<RANK,true> >(base);}


} // structure



#endif // STRUCTURE_STRUCTURE_H_INCLUDED
