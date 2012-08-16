// -*- C++ -*-
#ifndef ELEMENTS_COMPOSITE_SUBSYSTEM_H_INCLUDED
#define ELEMENTS_COMPOSITE_SUBSYSTEM_H_INCLUDED

#include "SubSystemFwd.h"

#include "Free.h"
#include "InteractionFwd.h"

#include "Structure.h"


namespace structure {


template<int RANK> 
class SubSystem
{
public:
  // const QuantumSystem<RANK>*const getQS() const {return qs_;}
  const Exact        <RANK>*const getEx() const {return ex_;} 
  const Hamiltonian  <RANK>*const getHa() const {return ha_;}
  const Liouvillean  <RANK>*const getLi() const {return li_;} 
  const Averaged     <RANK>*const getAv() const {return av_;}  

protected:
  SubSystem(const DynamicsBase* qs)
    : // qs_(qs),
      ex_(qse<RANK>(qs)),
      ha_(qsh<RANK>(qs)),
      li_(qsl<RANK>(qs)),
      av_(qsa<RANK>(qs))
  {}

  SubSystem() : ex_(), ha_(), li_(), av_() {}

private:
  // const QuantumSystem<RANK>* qs_;
  const Exact        <RANK>* ex_; 
  const Hamiltonian  <RANK>* ha_;
  const Liouvillean  <RANK>* li_; 
  const Averaged     <RANK>* av_;
  
};



template<int RANK> 
class SubSystemsInteraction : public SubSystem<RANK>
{
public:
  typedef class Interaction<RANK> Interaction;

  SubSystemsInteraction(const Interaction* ia) : SubSystem<RANK>(ia), ia_(ia) {}

  const Interaction*const get() const {return ia_;} 

private:
  const Interaction* ia_;

};



class SubSystemFree : public SubSystem<1>
{
public:
  SubSystemFree(const Free* free) : SubSystem<1>(free), free_(free) {}

  SubSystemFree() : free_() {}

  const Free*const get() const {return free_;}

private:
  const Free* free_;

};


} // structure


#endif // ELEMENTS_COMPOSITE_SUBSYSTEM_H_INCLUDED
