// -*- C++ -*-
#ifndef _STRUCTURE_H
#define _STRUCTURE_H

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


template<int RANK>
inline 
const Exact<RANK>*const 
qse(const QuantumSystem<RANK>* quantumSystem)
{return dynamic_cast<const Exact<RANK>*>(quantumSystem);}

template<int RANK>
inline 
const Hamiltonian<RANK,TWO_TIME>*const 
qsh(const QuantumSystem<RANK>* quantumSystem)
{return dynamic_cast<const Hamiltonian<RANK,TWO_TIME>*>(quantumSystem);}

template<int RANK>
inline 
const Liouvillean<RANK,true>*const 
qsl(const QuantumSystem<RANK>* quantumSystem)
{return dynamic_cast<const Liouvillean<RANK,true>*>(quantumSystem);}

template<int RANK>
inline 
const Averaged<RANK,true>*const 
qsa(const QuantumSystem<RANK>* quantumSystem)
{return dynamic_cast<const Averaged<RANK,true>*>(quantumSystem);}



template<int RANK>
inline 
const Exact<RANK>*const 
qse(const DynamicsBase* base)
{return dynamic_cast<const Exact<RANK>*>(base);}

template<int RANK>
inline 
const Hamiltonian<RANK,TWO_TIME>*const 
qsh(const DynamicsBase* base)
{return dynamic_cast<const Hamiltonian<RANK,TWO_TIME>*>(base);}

template<int RANK>
inline 
const Liouvillean<RANK,true>*const 
qsl(const DynamicsBase* base)
{return dynamic_cast<const Liouvillean<RANK,true>*>(base);}

template<int RANK>
inline 
const Averaged<RANK,true>*const 
qsa(const DynamicsBase* base)
{return dynamic_cast<const Averaged<RANK,true>*>(base);}


} // structure



#endif // _STRUCTURE_H
