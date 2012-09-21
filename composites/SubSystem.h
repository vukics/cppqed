// -*- C++ -*-
#ifndef ELEMENTS_COMPOSITES_SUBSYSTEM_H_INCLUDED
#define ELEMENTS_COMPOSITES_SUBSYSTEM_H_INCLUDED

#include "SubSystemFwd.h"

#include "Free.h"
#include "Interaction.h"

#include "Structure.h"


namespace composite {


template<int RANK> 
class SubSystemsInteraction : public structure::QuantumSystemWrapper<RANK>
{
public:
  typedef typename structure::Interaction<RANK>::Ptr InteractionPtr;

  SubSystemsInteraction(InteractionPtr ia) : structure::QuantumSystemWrapper<RANK>(ia), ia_(ia) {}

  const InteractionPtr get() const {return ia_;} 

private:
  InteractionPtr ia_;

};



class SubSystemFree : public structure::QuantumSystemWrapper<1>
{
public:
  typedef structure::Free::SmartPtr FreePtr;

  SubSystemFree(FreePtr free) : structure::QuantumSystemWrapper<1>(free,true), free_(free) {}

  SubSystemFree() : free_() {}

  const FreePtr get() const {return free_;}

private:
  FreePtr free_;

};


} // composite


#endif // ELEMENTS_COMPOSITES_SUBSYSTEM_H_INCLUDED
