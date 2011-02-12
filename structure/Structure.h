// -*- C++ -*-
#ifndef _STRUCTURE_H
#define _STRUCTURE_H

#include "QuantumSystemFwd.h"
#include "DynamicsBaseFwd.h"

#include "Exact.h"
#include "Hamiltonian.h"
#include "Liouvillean.h"
#include "Averaged.h"


namespace structure {


template<int RANK>
inline 
const Exact      <RANK>*const 
qse(const QuantumSystem<RANK>* quantumSystem)
{return dynamic_cast<const Exact      <RANK>*>(quantumSystem);}

template<int RANK>
inline 
const Hamiltonian<RANK>*const 
qsh(const QuantumSystem<RANK>* quantumSystem)
{return dynamic_cast<const Hamiltonian<RANK>*>(quantumSystem);}

template<int RANK>
inline 
const Liouvillean<RANK>*const 
qsl(const QuantumSystem<RANK>* quantumSystem)
{return dynamic_cast<const Liouvillean<RANK>*>(quantumSystem);}

template<int RANK>
inline 
const Averaged   <RANK>*const 
qsa(const QuantumSystem<RANK>* quantumSystem)
{return dynamic_cast<const Averaged   <RANK>*>(quantumSystem);}



template<int RANK>
inline 
const Exact      <RANK>*const 
qse(const DynamicsBase* base)
{return dynamic_cast<const Exact      <RANK>*>(base);}

template<int RANK>
inline 
const Hamiltonian<RANK>*const 
qsh(const DynamicsBase* base)
{return dynamic_cast<const Hamiltonian<RANK>*>(base);}

template<int RANK>
inline 
const Liouvillean<RANK>*const 
qsl(const DynamicsBase* base)
{return dynamic_cast<const Liouvillean<RANK>*>(base);}

template<int RANK>
inline 
const Averaged   <RANK>*const 
qsa(const DynamicsBase* base)
{return dynamic_cast<const Averaged   <RANK>*>(base);}


} // structure



#endif // _STRUCTURE_H
