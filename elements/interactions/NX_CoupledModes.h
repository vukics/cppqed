// -*- C++ -*-
#ifndef   NX_COUPLED_MODES_INCLUDED
#define   NX_COUPLED_MODES_INCLUDED

#include "NX_CoupledModesFwd.h"

#include "Mode_.h"

#include "Interaction.h"

namespace nxcoupledmodes {

class Base: public structure::Interaction<2>, public structure::TridiagonalHamiltonian<2,true>
{
public:
  Base(const ModeBase*, const ModeBase*, double u);

};

} // nxcoupledmodes



#define BIG_NAMESPACE_NAME             nxcoupledmodes
#define BIG_CLASS_NAME                 NX_CoupledModes
#define BIG_ADDITIONAL_PARAMETERS      , double u
#define BIG_ADDITIONAL_PARAMETERS_PASS ,u

#include "details/BinaryInteractionGenerator.h"

#endif // NX_COUPLED_MODES_INCLUDED
