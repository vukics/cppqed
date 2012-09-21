// -*- C++ -*-
#ifndef STRUCTURE_STRUCTURE_H_INCLUDED
#define STRUCTURE_STRUCTURE_H_INCLUDED

#include "StructureFwd.h"

#include "QuantumSystem.h"
#include "DynamicsBase.h"

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



template<int RANK> 
class QuantumSystemWrapper
{
public:
  typedef typename QuantumSystem<RANK>::Ptr QuantumSystemPtr;
  typedef typename Exact        <RANK>::Ptr ExactPtr;
  typedef typename Hamiltonian  <RANK>::Ptr HamiltonianPtr;
  typedef typename Liouvillean  <RANK>::Ptr LiouvilleanPtr;
  typedef typename Averaged     <RANK>::Ptr AveragedPtr;

  const QuantumSystemPtr getQS() const {return qs_;}
  const ExactPtr         getEx() const {return ex_;} 
  const HamiltonianPtr   getHa() const {return ha_;}
  const LiouvilleanPtr   getLi() const {return li_;} 
  const AveragedPtr      getAv() const {return av_;}  

  /*
  QuantumSystemPtr getQS() {return qs_;}
  ExactPtr getEx() {return ex_;} 
  HamiltonianPtr getHa() {return ha_;}
  LiouvilleanPtr getLi() {return li_;} 
  AveragedPtr getAv() {return av_;}  
  */

  QuantumSystemWrapper(DynamicsBase::Ptr qs)
    : qs_(dynamic_pointer_cast<const QuantumSystem<RANK> >(qs)),
      ex_(qse<RANK>(qs)),
      ha_(qsh<RANK>(qs)),
      li_(qsl<RANK>(qs)),
      av_(qsa<RANK>(qs))
  {}

  QuantumSystemWrapper(QuantumSystemPtr qs, bool noise=true)
    : qs_(qs),
      ex_(qse<RANK>(qs)),
      ha_(qsh<RANK>(qs)),
      li_(noise ? qsl<RANK>(qs) : LiouvilleanPtr()),
      av_(qsa<RANK>(qs))
  {}

protected:
  QuantumSystemWrapper() : qs_(), ex_(), ha_(), li_(), av_() {}

private:
  QuantumSystemPtr qs_;
  ExactPtr         ex_; 
  HamiltonianPtr   ha_;
  LiouvilleanPtr   li_; 
  AveragedPtr      av_;
  
};



} // structure



#endif // STRUCTURE_STRUCTURE_H_INCLUDED
