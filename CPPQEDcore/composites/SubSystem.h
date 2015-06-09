// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines composite::SubSystemFree and composite::SubSystemsInteraction}
#ifndef CPPQEDCORE_COMPOSITES_SUBSYSTEM_H_INCLUDED
#define CPPQEDCORE_COMPOSITES_SUBSYSTEM_H_INCLUDED

#include "SubSystemFwd.h"

#include "Free.h"
#include "Interaction.h"

#include "Structure.h"


namespace composite {


template<int RANK> 
class SubSystemsInteraction : public structure::QuantumSystemWrapper<RANK,false>
{
public:
  typedef typename structure::Interaction<RANK>::Ptr InteractionPtr;

  explicit SubSystemsInteraction(InteractionPtr ia) : structure::QuantumSystemWrapper<RANK,false>(ia), ia_(ia) {}

  const InteractionPtr get() const {return ia_;} 

private:
  InteractionPtr ia_;

};



class SubSystemFree : public structure::QuantumSystemWrapper<1,false>
{
public:
  typedef structure::Free::Ptr FreePtr;

  explicit SubSystemFree(FreePtr free) : structure::QuantumSystemWrapper<1,false>(free,true), free_(free) {}

  SubSystemFree() : free_() {}

  const FreePtr get() const {return free_;}

private:
  FreePtr free_;

};


} // composite


#endif // CPPQEDCORE_COMPOSITES_SUBSYSTEM_H_INCLUDED
