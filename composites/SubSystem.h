// -*- C++ -*-
#ifndef ELEMENTS_COMPOSITES_SUBSYSTEM_H_INCLUDED
#define ELEMENTS_COMPOSITES_SUBSYSTEM_H_INCLUDED

#include "SubSystemFwd.h"

#include "Free.h"
#include "Interaction.h"

#include "Structure.h"


namespace structure {


template<int RANK> 
class SubSystem
{
public:
  // const QuantumSystem<RANK>*const getQS() const {return qs_;}
  const typename Exact      <RANK>::Ptr getEx() const {return ex_;} 
  const typename Hamiltonian<RANK>::Ptr getHa() const {return ha_;}
  const typename Liouvillean<RANK>::Ptr getLi() const {return li_;} 
  const typename Averaged   <RANK>::Ptr getAv() const {return av_;}  

protected:
  SubSystem(DynamicsBase::Ptr qs)
    : // qs_(qs),
      ex_(qse<RANK>(qs)),
      ha_(qsh<RANK>(qs)),
      li_(qsl<RANK>(qs)),
      av_(qsa<RANK>(qs))
  {}

  SubSystem() : ex_(), ha_(), li_(), av_() {}

private:
  // const QuantumSystem<RANK>* qs_;
  typename Exact        <RANK>::Ptr ex_; 
  typename Hamiltonian  <RANK>::Ptr ha_;
  typename Liouvillean  <RANK>::Ptr li_; 
  typename Averaged     <RANK>::Ptr av_;
  
};



template<int RANK> 
class SubSystemsInteraction : public SubSystem<RANK>
{
public:
  typedef class Interaction<RANK> Interaction;

  SubSystemsInteraction(typename Interaction::Ptr ia) : SubSystem<RANK>(ia), ia_(ia) {}

  const typename Interaction::Ptr get() const {return ia_;} 

private:
  typename Interaction::Ptr ia_;

};



class SubSystemFree : public SubSystem<1>
{
public:
  SubSystemFree(Free::SmartPtr free) : SubSystem<1>(free), free_(free) {}

  SubSystemFree() : free_() {}

  const Free::SmartPtr get() const {return free_;}

private:
  Free::SmartPtr free_;

};


} // structure


#endif // ELEMENTS_COMPOSITES_SUBSYSTEM_H_INCLUDED
