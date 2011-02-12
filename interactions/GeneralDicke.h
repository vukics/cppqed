// -*- C++ -*-
#ifndef   _GENERAL_DICKE_INCLUDED
#define   _GENERAL_DICKE_INCLUDED

#include "GeneralDickeFwd.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode.h"
#include "Spin.h"


namespace generaldicke {

class Base : public structure::Interaction<2>, public structure::TridiagonalHamiltonian<2,true>
{
public:
  Base(const ModeBase*, const SpinBase*, dcomp u, dcomp y);

};

} // generaldicke


#define BIG_NAMESPACE_NAME             generaldicke
#define BIG_CLASS_NAME                 GeneralDicke
#define BIG_ADDITIONAL_PARAMETERS      , dcomp u, dcomp y
#define BIG_ADDITIONAL_PARAMETERS_PASS , u, y

#include "details/BinaryInteractionGenerator.h"




#endif // _GENERAL_DICKE_INCLUDED
