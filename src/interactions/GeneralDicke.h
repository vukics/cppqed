// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED

#include "GeneralDickeFwd.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode_.h"
#include "Spin.h"


namespace generaldicke {

class Base : public structure::Interaction<2>, public structure::TridiagonalHamiltonian<2,true>
{
public:
  Base(mode::Ptr, spin::Ptr, dcomp u, dcomp y);

};

} // generaldicke


#define BIG_NAMESPACE_NAME             generaldicke
#define BIG_CLASS_NAME                 GeneralDicke
#define BIG_ADDITIONAL_PARAMETERS      , dcomp u, dcomp y
#define BIG_ADDITIONAL_PARAMETERS_PASS , u, y

#include "details_BinaryInteractionGenerator.h"




#endif // ELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED
