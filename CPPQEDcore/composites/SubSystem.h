/// \briefFile{Defines composite::SubSystemFree and composite::SubSystemsInteraction}
#ifndef COMPOSITES_SUBSYSTEM_H_INCLUDED
#define COMPOSITES_SUBSYSTEM_H_INCLUDED

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


#endif // COMPOSITES_SUBSYSTEM_H_INCLUDED
