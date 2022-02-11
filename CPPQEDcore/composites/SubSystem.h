// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines composite::SubSystemFree and composite::SubSystemsInteraction}
#ifndef CPPQEDCORE_COMPOSITES_SUBSYSTEM_H_INCLUDED
#define CPPQEDCORE_COMPOSITES_SUBSYSTEM_H_INCLUDED

#include "Interaction.h"

#include "Structure.h"


namespace composite {


using ::structure::FreePtr, ::structure::InteractionPtr;

  
template<int RANK> 
class SubSystemsInteraction : public structure::QuantumSystemWrapper<RANK>
{
public:
  explicit SubSystemsInteraction(InteractionPtr<RANK> ia) : structure::QuantumSystemWrapper<RANK>(ia), ia_(ia) {}

  const InteractionPtr<RANK> get() const {return ia_;} 

private:
  InteractionPtr<RANK> ia_;

};



class SubSystemFree : public structure::QuantumSystemWrapper<1>
{
public:
  explicit SubSystemFree(FreePtr free) : structure::QuantumSystemWrapper<1>(free,true), free_(free) {}

  SubSystemFree() : free_() {}

  const FreePtr get() const {return free_;}

private:
  FreePtr free_;

};


} // composite


#endif // CPPQEDCORE_COMPOSITES_SUBSYSTEM_H_INCLUDED
